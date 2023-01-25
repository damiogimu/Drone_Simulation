#include "../my_header.h"

#define PLOT_FILE "DATA/result_line_path"
#define DESIRE_FILE "DATA/desire_path"

#define DATA_LABEL "angle [rad]"
#define X_MIN 0.0
#define X_ITV 1.0

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

//	fprintf(gp, "set \n");
	fprintf(gp, "set terminal qt font 'Times New Roman, 26'\n");
	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set xzeroaxis\n");
	fprintf(gp, "set yzeroaxis\n");
	fprintf(gp, "set xtics %f, %f, %f\n", X_MIN, X_ITV, TIME);
	fprintf(gp, "set xlabel 't [s]'\n");
	fprintf(gp, "set ylabel '%s'\n", DATA_LABEL);
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key spacing 1.5 font ',16'\n");

	fprintf(gp, "plot [0.0:%f][-1.5:1.5]\
			'%s' u 1:6 w l lw 2 lt 7 lc 'blue' title 'θ', \
			'%s' u 1:2 w l lt 0 lc 'grey' title 'desire x acc', \
			'%s' u 1:15 w l lt 0 lc 'red' title 'desire θ', \
			\n", TIME, PLOT_FILE, DESIRE_FILE, DESIRE_FILE);

	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
