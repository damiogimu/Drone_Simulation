#include "../my_header.h"

#define PLOT_FILE1 "DATA/desire_path"
#define PLOT_FILE2 "DATA/result_line_path"

#define X_MIN 0.0

#define Y_MIN -2.0
#define Y_MAX 4.0

int main(void)
{
	FILE *gp;
	char buf[5];

	if ((gp = popen("/usr/local/bin/gnuplot", "w")) == NULL)
	{
		fprintf(stderr, "file does not exist");
		return (-1);
	}

//	fprintf(gp, "set \n");
	fprintf(gp, "set terminal qt font 'Times New Roman, 26'\n");
	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set grid\n");
	fprintf(gp, "set xtics %f, 1, %f\n", X_MIN, TIME);
	fprintf(gp, "set ytics %f, 1, %f\n", Y_MIN, Y_MAX);
	fprintf(gp, "set xlabel 't [s]'\n");
	fprintf(gp, "set ylabel 'Acceleration [m/s^2]'\n");
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key spacing 1.5 font ',16'\n");

	fprintf(gp, "plot [%f:%f] [%f:%f] '%s' u 1:8  w l lt 0 lc 'red' title 'x velo desire'\n", X_MIN, TIME, Y_MIN, Y_MAX, PLOT_FILE1);
	fprintf(gp, "replot '%s' u 1:10  w l lc 'blue' title 'x velo real'\n", PLOT_FILE2);

	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}