#include "my_header.h"

#define PLOT_FILE "DATA/result_lag" 

#define X_MIN 0.0
#define X_MAX 5.0

#define Y_MIN 0.0
#define Y_MAX 3.0

#define REF_DATA 8 

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
	fprintf(gp, "set terminal qt font 'Times New Roman, 20'\n");
	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set xzeroaxis\n");
	fprintf(gp, "set xtics \n");
	fprintf(gp, "set ytics \n");
	fprintf(gp, "set xlabel 'T [s]'\n");
	fprintf(gp, "set ylabel 'z [m]'\n");
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key font ',16'\n");

	fprintf(gp, "plot [0.0:%f][-1.0:1.0] '%s' u 1:%d w l notitle\n", TIME, PLOT_FILE, REF_DATA);
	
	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
