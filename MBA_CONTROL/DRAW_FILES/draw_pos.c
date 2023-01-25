#include "../my_header.h"

#define PLOT_FILE1 "DATA/result_line_path"
#define PLOT_FILE2 "DATA/desire_path"

#define X_MIN 0.0
#define X_ITV 5.0

#define Y_MIN -1.0
#define Y_MAX 3.0
#define Y_ITV 1.0

#define X_FLAG 1

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
	fprintf(gp, "set xtics %f, %f, %f\n", X_MIN, X_ITV, TIME);
	fprintf(gp, "set ytics %f, %f, %f\n", Y_MIN, Y_ITV, Y_MAX);
	fprintf(gp, "set xlabel 't [s]'\n");
	fprintf(gp, "set ylabel 'Position [m]'\n");
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key spacing 1.5 font ',16'\n");

	if (X_FLAG == 1)
	{
		fprintf(gp, "plot [%f:%f] [%f:%f] \
			'%s' u 1:2 w l lw 2.5 lt 7 lc 'blue' title 'x pos', \
			'%s' u 1:11 w l lt 0 lc 'red' title 'desire', \
			\n", X_MIN, TIME, Y_MIN, Y_MAX, PLOT_FILE1, PLOT_FILE2);
	}
	else
	{
		fprintf(gp, "plot [%f:%f] [%f:%f] \
			'%s' u 1:3 w l lw 2.5 lt 7 lc 'green' title 'y pos', \
			'%s' u 1:12 w l lt 0 lc 'red' title 'desire', \
			\n", X_MIN, TIME, Y_MIN, Y_MAX, PLOT_FILE1, PLOT_FILE2);
	}

	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
