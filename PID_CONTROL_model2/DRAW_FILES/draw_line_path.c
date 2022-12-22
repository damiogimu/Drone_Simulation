#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define PLOT_FILE "DATA/result_line_path"

#define X_MIN 0.0
#define X_MAX 4.0

#define Y_MIN 0.0
#define Y_MAX 4.0

#define Z_MIN 1.2
#define Z_MAX 3.2

int main(void)
{
	FILE *gp;
	char buf[5];

	if ((gp = popen("/usr/local/bin/gnuplot", "w")) == NULL)
	{
		fprintf(stderr, "file does not exist");
		return (-1);
	}

	if (X_MIN != X_MAX)
		fprintf(gp, "set xrange[%f:%f]\n", X_MIN, X_MAX);
	if (Y_MIN != Y_MAX)
		fprintf(gp, "set yrange[%f:%f]\n", Y_MIN, Y_MAX);
	if (Z_MIN != Z_MAX)
		fprintf(gp, "set zrange[%f:%f]\n", Z_MIN, Z_MAX);

//	fprintf(gp, "set \n");
	fprintf(gp, "set terminal qt font 'Times New Roman, 20'\n");
	fprintf(gp, "set view 72.0, 347.0\n");
	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set grid\n");
	fprintf(gp, "set ticslevel 0\n");
//	fprintf(gp, "set xzeroaxis\n");
//	fprintf(gp, "set yzeroaxis\n");
//	fprintf(gp, "set zzeroaxis\n");
	fprintf(gp, "set xtics %f, 1, %f\n", X_MIN, X_MAX);
	fprintf(gp, "set ytics %f, 1, %f\n", Y_MIN, Y_MAX-1.0);
	fprintf(gp, "set ztics %f, 0.2, %f\n", Z_MIN, Z_MAX);
	fprintf(gp, "set xlabel 'x [m]'\n");
	fprintf(gp, "set ylabel 'y [m]'\n");
	fprintf(gp, "set zlabel 'z [m]'\n");
	fprintf(gp, "set border lw 1\n");
	fprintf(gp, "set key spacing 1.5 font ',16'\n");

//	fprintf(gp, "set terminal png\n");
//	fprintf(gp, "set output 'IMAGES/line_path.png'\n");	

	fprintf(gp, "splot '%s' u 2:3:4 w l lw 1 lt 7 lc 'blue' notitle\n", PLOT_FILE);
	fprintf(gp, "replot '%s' u 18:19:20 w l lw 0.5 lt 7 lc 'green' notitle\n", PLOT_FILE);
	fprintf(gp, "replot '%s' u 2:3:4 w l lw 0.8 lt 7 dt (10,15) lc 'red' notitle\n", "DATA/desire_line_path");

	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
