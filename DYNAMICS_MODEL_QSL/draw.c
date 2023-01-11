#include "my_header.h"

#ifndef LAG_F
# define LAG_F 0
#endif
#ifndef ERROR_F
# define ERROR_F 0
#endif
#ifndef ACC_F
# define ACC_F 0
#endif

#define NEWTON_FILE "./RESULT/result_newton"
#define LAG_FILE "./RESULT/result_lag"
#define ERROR_FILE "./RESULT/result_error"

#define X_MIN 0.0
#define X_MAX 5.0

#define Y_MIN 0.0
#define Y_MAX 3.0

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
	fprintf(gp, "set xtics \n");
	fprintf(gp, "set ytics \n");
	fprintf(gp, "set xlabel 'T [s]'\n");
	fprintf(gp, "set ylabel 'z [m]'\n");
	if (ERROR_F == 1)
		fprintf(gp, "set ylabel 'Error [m]'\n");
	else if (ACC_F == 0)
		fprintf(gp, "set ylabel 'Position [m]'\n");
	else
		fprintf(gp, "set ylabel 'ACC [m/s^2]'\n");
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key font ',16'\n");


	if (ACC_F == 1)
	{
			fprintf(gp, "plot [0:%f][0:0.4] \
					'%s' u 1:18 w l notitle, \
					'%s' u 1:19 w l notitle, \
					\n", TIME, LAG_FILE, LAG_FILE);
	}
	else
	{
		if (ERROR_F == 1)
		{
			fprintf(gp, "plot [0:%f][-1:1] \
					'%s' u 1:4 w l notitle, \
					\n", TIME, ERROR_FILE);
		}
		else
		{
			if (LAG_F == 0)
			{
				fprintf(gp, "plot [0:%f][0:4] \
						'%s' u 1:4 w l notitle, \
						\n", TIME, NEWTON_FILE);
			}	
			else
			{
				fprintf(gp, "plot [0:%f][0:4] \
						'%s' u 1:4 w l notitle, \
						\n", TIME, LAG_FILE);
			}
		}
	}
	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
