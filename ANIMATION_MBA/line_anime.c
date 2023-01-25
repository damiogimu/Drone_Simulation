#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "../MBA_CONTROL/my_header.h"

#define INTERVAL 0.05

#define X_INIT_MIN -1.0
#define X_INIT_MAX 3.0
#define X_ITV 1.0

#define Y_INIT_MIN -1.0
#define Y_INIT_MAX 3.0
#define Y_ITV 1.0

#define Z_INIT_MIN 0.0
#define Z_INIT_MAX 3.0
#define Z_ITV 1.0

#ifndef LOAD_F
# define LOAD_F 0
#endif

int main(void)
{
	FILE *gp;
	int i = 0;
	double t = 0.0;

	if ((gp = popen("/usr/local/bin/gnuplot", "w")) == NULL)
	{
		fprintf(stderr, "gnuplot does not exist\n");
		return (-1);
	}
	
//	fprintf(gp, "\n");
	fprintf(gp, "reset\n");
	fprintf(gp, "set nokey\n");
	fprintf(gp, "set grid\n");
	fprintf(gp, "set ticslevel 0\n");
	fprintf(gp, "set xrange [%f:%f]\n", X_INIT_MIN, X_INIT_MAX);
	fprintf(gp, "set yrange [%f:%f]\n", Y_INIT_MIN, Y_INIT_MAX);
	fprintf(gp, "set zrange [%f:%f]\n", Z_INIT_MIN, Z_INIT_MAX);
	fprintf(gp, "set xtics %f, %f, %f\n", X_INIT_MIN, X_ITV, X_INIT_MAX);
	fprintf(gp, "set ytics %f, %f, %f\n", Y_INIT_MIN, Y_ITV, Y_INIT_MAX);
	fprintf(gp, "set ztics %f, %f, %f\n", Z_INIT_MIN, Z_ITV, Z_INIT_MAX);
	fprintf(gp, "set xlabel 'x [m]' offset 0.0,-1.5,0.0\n");
	fprintf(gp, "set ylabel 'y [m]' offset -1.0,-1.0,0.0\n");
	fprintf(gp, "set zlabel 'z [m]' offset 1.5,0.0,0.0\n");
	fprintf(gp, "set border lw 1\n");
	fprintf(gp, "set view 60,345\n");
//	fprintf(gp, "set view 90,360\n");

	fprintf(gp, "set term gif animate delay 5 font 'Times New Roman, 20'\n");
	fprintf(gp, "set output 'GIF/MBA_PATH.gif'\n");

	while (t <= TIME)
	{
		if (LOAD_F == 1)
		{	
			fprintf(gp, "splot \
					'DATA/xrotor_DATA' index %d u 2:3:4 w l lc 'blue', \
					'DATA/yrotor_DATA' index %d u 2:3:4 w l lc 'blue', \
					'DATA/cable_DATA' index %d u 2:3:4 w l lc 'black', \
					'DATA/path_DATA' every ::0::%d u 2:3:4 w l lw 0.7 lc 'blue',\
					'DATA/path_DATA' every ::%d::%d u 2:3:4 with points pt 13 ps 1.0 lc 'blue',\
					'DATA/path_DATA' every ::0::%d u 8:9:10 w l lw 0.7 lc 'green',\
					'DATA/path_DATA' every ::%d::%d u 8:9:10 with points pt 7 ps 1.5 lc 'green',\
					'DATA/path_DATA' u 11:12:13 w l lt 0 lc 'red',\
					\n", i, i, i, i, i, i, i, i, i);
		}
		else
		{
			fprintf(gp, "splot \
					'DATA/xrotor_DATA' index %d u 2:3:4 w l lc 'blue', \
					'DATA/yrotor_DATA' index %d u 2:3:4 w l lc 'blue', \
					'DATA/path_DATA' every ::0::%d u 2:3:4 w l lw 0.7 lc 'blue',\
					'DATA/path_DATA' every ::%d::%d u 2:3:4 with points pt 13 ps 1.0 lc 'blue',\
					'DATA/path_DATA' u 11:12:13 w l lt 0 lc 'red',\
					\n", i, i, i, i, i);
		}
		i++;
		t += INTERVAL;
	}

	fflush(gp);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
