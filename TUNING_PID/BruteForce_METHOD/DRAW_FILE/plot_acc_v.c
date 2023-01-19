#include "../my_header.h"

#define X_MIN 0.0

#define Y_MIN -6.0
#define Y_MAX 6.0

int main(void)
{
	FILE *gp;
	char *filename;
	int Kpt, Kdt, Kit, Kpr, Kdr, Kir;
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

	Kpt = MIN_Pt_GAIN;
	while (Kpt <= MAX_Pt_GAIN)
	{
		Kdt = MIN_Dt_GAIN;
		while (Kdt <= MAX_Dt_GAIN)
		{
			Kit = MIN_It_GAIN;
			while (Kit <= MAX_It_GAIN)
			{
				Kpr = MIN_Pr_GAIN;
				while (Kpr <= MAX_Pr_GAIN)
				{
					Kdr = MIN_Dr_GAIN;
					while (Kdr <= MAX_Dr_GAIN)
					{
						Kir = MIN_Ir_GAIN;
						while (Kir <= MAX_Ir_GAIN)
						{
							// --- Judge file exist --- //
							sprintf(filename, "DATA/gain%d_%d_%d_%d_%d_%d", Kpt, Kdt, Kit, Kpr, Kdr, Kir);
							if ((fd = open(filename, O_RDONLY)) != -1)
							{
								fprintf(gp, "plot \
										'./DATA/gain%d_%d_%d_%d_%d_%d' u 1:2 w l lc 'grey' title '%d %d %d %d %d %d x_acc',\
										'./DATA/gain%d_%d_%d_%d_%d_%d' u 1:4 w l lt 0 lc 'orange' title 'y_acc',\
										\n", Kpt, Kdt, Kit, Kpr, Kdr, Kir, Kpt, Kdt, Kit, Kpr, Kdr, Kir, Kpt, Kdt, Kit, Kpr, Kdr, Kir);
							}
							// ------------------------ //
							close(fd);
							Kir += Ir_ITV;
						}
						Kdr += Dr_ITV;
					}
					Kpr += Pr_ITV;
				}
				Kit += It_ITV;
			}
			Kdt += Dt_ITV;
		}
		Kpt += Pt_ITV;
	}

	fflush(gp);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
