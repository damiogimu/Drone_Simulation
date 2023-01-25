#include "../my_header.h"

#define X_MIN 0.0

#define Y_MIN -3.0
#define Y_MAX 3.0

int main(void)
{
	FILE *gp;
	char *filename;
	int fd;
	int Kpt, Kdt, Kpr, Kdr, Kir;

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
	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set xtics %f, 0.5, %f\n", X_MIN, TIME);
	fprintf(gp, "set xlabel 't [s]'\n");
	fprintf(gp, "set ylabel 'Acceleration [m/s^2]'\n");
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set term gif animate delay 10 font 'Times New Roman, 20'\n");
	fprintf(gp, "set output 'PLOT_DATA/EACH_GAIN_ACC_P.gif'\n");

	Kpt = MIN_Kpt_GAIN;
	while (Kpt <= MAX_Kpt_GAIN)
	{
		Kdt = MIN_Kdt_GAIN;
		while (Kdt <= MAX_Kdt_GAIN)
		{
			Kpr = MIN_P_GAIN;
			while (Kpr <= MAX_P_GAIN)
			{
				Kdr = MIN_D_GAIN;
				while (Kdr <= MAX_D_GAIN)
				{
					Kir = MIN_I_GAIN;
					while (Kir <= MAX_I_GAIN)
					{
						sprintf(filename, "DATA/gain%d_%d_%d_%d_%d", Kpt, Kdt, Kpr, Kdr, Kir);
						if ((fd = open(filename, O_RDONLY)) != -1)
						{
							fprintf(gp, "plot \
									'./DATA/gain%d_%d_%d_%d_%d' u 1:2 w l lc 'grey' title '%d %d %d %d %d φ',\
									'./DATA/gain%d_%d_%d_%d_%d' u 1:4 w l lt 0 lc 'orange' title 'θ',\
									\n", Kpt, Kdt, Kpr, Kdr, Kir, Kpt, Kdt, Kpr, Kdr, Kir, Kpt, Kdt, Kpr, Kdr, Kir);
						}
						close(fd);
						Kir += I_ITV;
					}
					Kdr += D_ITV;
				}
				Kpr += P_ITV;
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
