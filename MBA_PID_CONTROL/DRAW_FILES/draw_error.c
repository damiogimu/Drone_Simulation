#include "../my_header.h"

#define PLOT_FILE1 "DATA/result_line_error"

#define X_MIN 0.0

#define Y_MIN -1.0
#define Y_MAX 1.0

int main(void)
{
	FILE *gp;
	char buf[5];

	if ((gp = popen("/usr/local/bin/gnuplot", "w")) == NULL)
	{
		fprintf(stderr, "file does not exist");
		return (-1);
	}

	if (X_MIN != TIME)
		fprintf(gp, "set xrange[%f:%f]\n", X_MIN, TIME);
	if (Y_MIN != Y_MAX)
		fprintf(gp, "set yrange[%f:%f]\n", Y_MIN, Y_MAX);

//	fprintf(gp, "set \n");
	fprintf(gp, "set terminal qt font 'Times New Roman, 26'\n");
	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set grid\n");
	fprintf(gp, "set xzeroaxis\n");
	fprintf(gp, "set yzeroaxis\n");
	fprintf(gp, "set xtics %f, 1, %f\n", X_MIN, TIME);
//	fprintf(gp, "set ytics %f, 100, %f\n", Y_MIN, Y_MAX);
	fprintf(gp, "set xlabel font 'Times New Roman, 35' 't [s]'\n");
	fprintf(gp, "set ylabel font 'Times New Roman, 35' 'Acceleration [m/s^2]'\n");
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key spacing 1.0\n");

	fprintf(gp, "plot '%s' u 1:2 w l lw 2 lt 7 lc 'red' title 'x'\n", PLOT_FILE1);
	fprintf(gp, "replot '%s' u 1:3 w l lw 1 lt 7 lc 'green' title 'y'\n", PLOT_FILE1);
	fprintf(gp, "replot '%s' u 1:4 w l lw 2 lt 7 lc 'blue' title 'z'\n", PLOT_FILE1);

	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
