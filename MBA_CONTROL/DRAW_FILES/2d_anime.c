#include "../my_header.h"

#define INTERVAL 0.05

#define X_INIT_MIN -1.0
#define X_INIT_MAX 3.0
#define X_ITV 1.0

#define Y_INIT_MIN 0.0
#define Y_INIT_MAX 3.0
#define Y_ITV 0.5

int main(void)
{
	FILE *gp;
	int i = 0;
	double t;

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
	fprintf(gp, "set xtics %f, %f, %f\n", X_INIT_MIN, X_ITV, X_INIT_MAX);
	fprintf(gp, "set ytics %f, %f, %f\n", Y_INIT_MIN, Y_ITV, Y_INIT_MAX - Y_ITV);
	fprintf(gp, "set xlabel 'x [m]'\n");
	fprintf(gp, "set ylabel 'z [m]'\n");
	fprintf(gp, "set border lw 1\n");
	fprintf(gp, "set view 90,360\n");

	fprintf(gp, "set term gif animate delay 5 font 'Times New Roman, 20'\n");
	fprintf(gp, "set output 'GIF/2D_ANIMATION.gif'\n");

	t = 0.0;
	while (t <= TIME)
	{
		fprintf(gp, "plot \
				'DATA/xrotor_DATA' index %d u 2:4 w l lw 2.0 lc 'blue', \
				'DATA/cable_DATA' index %d u 2:4 w l lc 'black', \
				'DATA/path_DATA' every ::%d::%d u 2:4 with points pt 13 ps 1.5 lc 'blue',\
				'DATA/path_DATA' every ::%d::%d u 8:10 with points pt 7 ps 2.0 lc 'green',\
				'DATA/path_DATA' u 11:13 w l lw 1.5 lt 0 lc 'red' \
				\n", i, i, i, i, i, i);
/*
		fprintf(gp, "plot \
				'DATA/xrotor_DATA' index %d u 2:4 w l lc 'blue', \
				'DATA/cable_DATA' index %d u 2:4 w l lc 'black', \
				'DATA/path_DATA' every ::0::%d u 2:4 w l lw 0.7 lc 'blue',\
				'DATA/path_DATA' every ::%d::%d u 2:4 with points pt 13 ps 1.0 lc 'blue',\
				'DATA/path_DATA' every ::0::%d u 8:10 w l lw 0.7 lc 'green',\
				'DATA/path_DATA' every ::%d::%d u 8:10 with points pt 7 ps 1.5 lc 'green',\
				'DATA/path_DATA' u 11:13 w l lt 0 lc 'red',\
				\n", i, i, i, i, i, i, i, i);
*/
		i++;
		t += INTERVAL;
	}

	fflush(gp);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}