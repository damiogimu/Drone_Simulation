#include "../my_header.h"

#define INTERVAL 0.05

#define X_INIT_MIN -0.5
#define X_INIT_MAX 3.0
#define X_ITV 0.5

#define Y_INIT_MIN -0.5
#define Y_INIT_MAX 3.0
#define Y_ITV 0.5

#define Z_INIT_MIN 0.0
#define Z_INIT_MAX 2.5
#define Z_ITV 0.5

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
//	fprintf(gp, "set xtics %f, %f, %f\n", X_INIT_MIN, X_ITV, X_INIT_MAX);
//	fprintf(gp, "set ytics %f, %f, %f\n", Y_INIT_MIN, Y_ITV, Y_INIT_MAX - Y_ITV);
	fprintf(gp, "set xtics (0,1.0,2.0,3.0) offset 0.0,-0.5,0.0\n");
	fprintf(gp, "set ytics (0,1.0,2.0,3.0) offset -1.0,-0.5,0.0\n");
	fprintf(gp, "set ztics %f, %f, %f\n", Z_INIT_MIN, Z_ITV, Z_INIT_MAX);
	fprintf(gp, "set xlabel 'x [m]' offset 0.0,-1.5,0.0\n");
	fprintf(gp, "set ylabel 'y [m]' offset -2.0,-2.0,0.0\n");
	fprintf(gp, "set zlabel 'z [m]' offset 1.5,0.0,0.0\n");
	fprintf(gp, "set border lw 1\n");
	fprintf(gp, "set view 70,345\n");

	fprintf(gp, "set term gif animate delay 5 font 'Times New Roman, 20'\n");
	fprintf(gp, "set output 'GIF/3d_GIF.gif'\n");

	while (t <= TIME)
	{
		fprintf(gp, "splot \
				'DATA/GIF_xrotor' index %d u 2:3:4 w l lw 2.0 lc 'blue', \
				'DATA/GIF_yrotor' index %d u 2:3:4 w l lw 2.0 lc 'blue', \
				'DATA/GIF_cable' index %d u 2:3:4 w l lc 'black', \
				'DATA/GIF_path' every ::%d::%d u 2:3:4 with points pt 13 ps 1.5 lc 'blue',\
				'DATA/GIF_path' every ::%d::%d u 8:9:10 with points pt 7 ps 2.0 lc 'green',\
				'DATA/GIF_path' u 11:12:13 w l lw 1.5 lt 0 lc 'red',\
				\n", i, i, i, i, i, i, i);
/*
		fprintf(gp, "splot \
				'DATA/GIF_xrotor' index %d u 2:3:4 w l lc 'blue', \
				'DATA/GIF_yrotor' index %d u 2:3:4 w l lc 'blue', \
				'DATA/GIF_cable' index %d u 2:3:4 w l lc 'black', \
				'DATA/GIF_path' every ::0::%d u 2:3:4 w l lw 0.7 lc 'blue',\
				'DATA/GIF_path' every ::%d::%d u 2:3:4 with points pt 13 ps 1.0 lc 'blue',\
				'DATA/GIF_path' every ::0::%d u 8:9:10 w l lw 0.7 lc 'green',\
				'DATA/GIF_path' every ::%d::%d u 8:9:10 with points pt 7 ps 1.5 lc 'green',\
				'DATA/GIF_path' u 11:12:13 w l lt 0 lc 'red',\
				\n", i, i, i, i, i, i, i, i, i);
*/
		i++;
		t += INTERVAL;
	}

	fflush(gp);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
