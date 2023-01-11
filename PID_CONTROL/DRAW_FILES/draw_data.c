#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "../my_val_define.h"

#define PLOT_FILE "DATA/result_line_path"

#ifndef REF_DATA
# define REF_DATA 8
#endif
#define DATA_LABEL "angle [rad]"

#define X_MIN 0.0
#define X_ITV 1.0

#define Y_MIN 0.0
#define Y_MAX 0.0

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
	fprintf(gp, "set xzeroaxis\n");
	fprintf(gp, "set yzeroaxis\n");
	fprintf(gp, "set xtics %f, %f, %f\n", X_MIN, X_ITV, TIME);
//	fprintf(gp, "set ytics %f, %f, %f\n", Y_MIN, Y_ITV, Y_MAX);
	fprintf(gp, "set xlabel 't [s]'\n");
	fprintf(gp, "set ylabel '%s'\n", DATA_LABEL);
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key spacing 1.5 font ',16'\n");

	fprintf(gp, "plot '%s' u 1:%d w l lw 1 lt 7 lc 'dark-orange' notitle\n", PLOT_FILE, REF_DATA);

	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
