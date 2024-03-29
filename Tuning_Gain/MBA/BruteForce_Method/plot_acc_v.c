#include "my_header.h"

#define X_MIN 0.0

#define Y_MIN -6.0
#define Y_MAX 6.0

int main(void)
{
	FILE *gp;
	char *filename;
	int Kpt, Kdt, Kpr, Kdr, Kir;
	int fd;

	filename = malloc(sizeof(char) * 80);

	if ((gp = popen("/usr/local/bin/gnuplot", "w")) == NULL)
	{
		fprintf(stderr, "file does not exist");
		return (-1);
	}

	if (X_MIN != TIME)
		fprintf(gp, "set xrange [%f:%f]\n", X_MIN, TIME);
	if (Y_MIN != Y_MAX)
		fprintf(gp, "set yrange [%f:%f]\n", Y_MIN, Y_MAX);

//	fprintf(gp, "set \n");
	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set xtics %f, 0.5, %f\n", X_MIN, TIME);
	fprintf(gp, "set xlabel 't [s]'\n");
	fprintf(gp, "set ylabel 'Acceleration [m/s^2]'\n");	
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set term gif animate delay 10 font 'Times New Roman, 20'\n");
	fprintf(gp, "set output 'PLOT_DATA/EACH_GAIN_ACC_V.gif'\n");

	Kpt = MIN_Kpt_GAIN;
	while (Kpt <= MAX_Kpt_GAIN)
	{
		Kdt = MIN_Kdt_GAIN;
		while (Kdt <= MAX_Kdt_GAIN)
		{
			Kpr = MIN_Kpr_GAIN;
			while (Kpr <= MAX_Kpr_GAIN)
			{
				Kdr = MIN_Kdr_GAIN;
				while (Kdr <= MAX_Kdr_GAIN)
				{
					// --- Judge file exist --- //
					sprintf(filename, "DATA/gain%d_%d_%d_%d", Kpt, Kdt, Kpr, Kdr);
					if ((fd = open(filename, O_RDONLY)) != -1)
					{
						fprintf(gp, "plot \
								'./DATA/gain%d_%d_%d_%d' u 1:2 w l lc 'grey' title '%d %d %d %d x_acc',\
								'./DATA/gain%d_%d_%d_%d' u 1:4 w l lt 0 lc 'orange' title 'y_acc',\
								\n", Kpt, Kdt, Kpr, Kdr, Kpt, Kdt, Kpr, Kdr, Kpt, Kdt, Kpr, Kdr);
					}
					// ------------------------ //
					close(fd);
					Kdr += Kdr_ITV;
				}
				Kpr += Kpr_ITV;
			}
			Kdt += Kdt_ITV;
		}
		Kpt += Kpt_ITV;
	}

	fflush(gp);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
