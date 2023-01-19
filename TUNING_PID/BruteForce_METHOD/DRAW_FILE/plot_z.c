#include "../my_header.h"

#define X_MIN 0.0

#define Y_MIN -6.0
#define Y_MAX 6.0

int main(void)
{
	FILE *gp;
	char *filename;
	int Kpt, Kdt, Kit;
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
	fprintf(gp, "set output 'PLOT_DATA/EACH_GAIN_Z.gif'\n");

	Kpt = MIN_Pz_GAIN;
	while (Kpt <= MAX_Pz_GAIN)
	{
		Kdt = MIN_Dz_GAIN;
		while (Kdt <= MAX_Dz_GAIN)
		{
			Kit = MIN_Iz_GAIN;
			while (Kit <= MAX_Iz_GAIN)
			{
				// --- Judge file exist --- //
				sprintf(filename, "Z_DATA/gain_%d_%d_%d", Kpt, Kdt, Kit);
				if ((fd = open(filename, O_RDONLY)) != -1)
				{
					fprintf(gp, "plot \
							'./Z_DATA/gain_%d_%d_%d' u 1:2 w l lc 'blue' title '%d %d %d x_acc',\
							'./Z_DATA/gain_%d_%d_%d' u 1:3 w l lt 0 lc 'red' title 'y_acc',\
							\n", Kpt, Kdt, Kit, Kpt, Kdt, Kit, Kpt, Kdt, Kit);
				}
				// ------------------------ //
				close(fd);
				Kit += Iz_ITV;
			}
			Kdt += Dz_ITV;
		}
		Kpt += Pz_ITV;
	}

	fflush(gp);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
