#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "../PID_CONTROL_model2/my_val_define.h"

#define INTERVAL 0.05
#define S_T 0.0

#define X_INIT_MIN -1.0
#define X_INIT_MAX 3.0
#define X_ITV 1.0

#define Y_INIT_MIN -1.0
#define Y_INIT_MAX 3.0
#define Y_ITV 1.0

#define Z_INIT_MIN -1.0
#define Z_INIT_MAX 2.0
#define Z_ITV 0.5

#ifndef LOAD_F
# define LOAD_F 0
#endif

int main(void)
{
	FILE *gp;
	int i = 0;
	double t = S_T;
	double x_min=X_INIT_MIN, x_max=X_INIT_MAX;
	double y_min=Y_INIT_MIN, y_max=Y_INIT_MAX;
	double z_min=Z_INIT_MIN, z_max=Z_INIT_MAX;
	int sc = S_T*20;

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
	fprintf(gp, "set ytics %f, %f, %f\n", Y_INIT_MIN, Y_ITV, Y_INIT_MAX-1);
	fprintf(gp, "set ztics %f, %f, %f\n", Z_INIT_MIN, Z_ITV, Z_INIT_MAX);
	fprintf(gp, "set xlabel 'x [m]' offset 0.0,-1.5,0.0\n");
	fprintf(gp, "set ylabel 'y [m]' offset -2.0,-2.0,0.0\n");
	fprintf(gp, "set zlabel 'z [m]' offset 1.5,0.0,0.0\n");
	fprintf(gp, "set border lw 1\n");
	fprintf(gp, "set view 45,345\n");

	fprintf(gp, "set term gif animate delay 5 font 'Times New Roman, 20'\n");
	fprintf(gp, "set output 'GIF/PID_PATH.gif'\n");

	while (t <= TIME)
	{
// ---- 途中 (ドローンの動きに対して軸の目盛をずらす) ---- //
/*
		fprintf(gp, "set xrange [%f:%f]\n", X_MIN, X_MAX);
		fprintf(gp, "set yrange [%f:%f]\n", Y_MIN, Y_MAX);
		fprintf(gp, "set zrange [%f:%f]\n", Z_MIN, Z_MAX);
		fprintf(gp, "set xtics %f, %d, %f offset 0.0,-0.5,0.0\n", X_MIN, (int)((fabs(X_MIN)+fabs(X_MAX))/X_AXIS_DIV), X_MAX);
		fprintf(gp, "set ytics %f, %d, %f\n", Y_MIN, (int)((fabs(Y_MIN)+fabs(Y_MAX))/Y_AXIS_DIV), Y_MAX);
		fprintf(gp, "set ztics %f+0.5*i*0.05*i*0.05, %d, %.1f+0.5*i*0.05*i*0.05\n", Z_MIN, (int)((fabs(Z_MIN)+fabs(Z_MAX))/Z_AXIS_DIV), Z_MAX);
*/
// --------------------------------------------------------- //
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
		sc++;
	}

	fflush(gp);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
