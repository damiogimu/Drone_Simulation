#include "../my_header.h"

#define PLOT_FILE "DATA/result_line_path"

#ifndef ATTI_F
# define ATTI_F 0
#endif
#ifndef REF_DATA
# define REF_DATA 9
#endif

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
	fprintf(gp, "set terminal qt font 'Arial, 20'\n");
	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set xzeroaxis\n");
	fprintf(gp, "set yzeroaxis\n");
	fprintf(gp, "set xtics %f, %f, %f\n", X_MIN, X_ITV, TIME);
	fprintf(gp, "set xlabel 't [s]'\n");
	fprintf(gp, "set ylabel '%s'\n", DATA_LABEL);
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key spacing 1.5 font ',16'\n");

	if (ATTI_F == 1)
		fprintf(gp, "plot [0.0:%f][-0.5:0.5]\
				'%s' u 1:5 w l lw 2 lt 7 lc 'blue' title 'φ', \
				'%s' u 1:6 w l lw 1	lt 7 lc 'light-green' title 'θ', \
				\n", TIME, PLOT_FILE, PLOT_FILE);
	else
		fprintf(gp, "plot [0.0:%f][-0.5:0.5]\
				'%s' u 1:8 w l lw 2 lt 7 lc 'dark-orange' title 'γ', \
				'%s' u 1:9 w l lw 1 lt 7 lc 'grey' title 'β', \
				\n", TIME, PLOT_FILE, PLOT_FILE);

	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
