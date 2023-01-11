#include "my_header.h"

#ifndef POS_F
# define POS_F 0
#endif

#ifndef LAG_F
# define LAG_F 0
#endif

#ifndef ERROR_F
# define ERROR_F 0
#endif

#define NEWTON_FILE "./RESULT/result_newton"
#define LAG_FILE "./RESULT/result_lag"
#define ERROR_FILE "./RESULT/result_error"

#define X_MIN 0.0
#define X_MAX 10.0

#define Y_MIN 0.0
#define Y_MAX 7.0

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
	if (POS_F == 1)
	{	
		if (ERROR_F == 1)
			fprintf(gp, "set ylabel 'Error [m]'\n");
	//	else
	//		fprintf(gp, "set ylabel 'Position [m]'\n");
	}
	else
	{
		if (ERROR_F == 1)
			fprintf(gp, "set ylabel 'Error [m]'\n");
		else
			fprintf(gp, "set ylabel 'angle [rad]'\n");
	}
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key font ',16'\n");

	if (ERROR_F == 1)
	{
		if (POS_F == 1)
		{
			fprintf(gp, "plot [0:10][-1:1]\
					'%s' u 1:4 w l lw 1.5  notitle, \
					\n", ERROR_FILE);
		}
		else
		{
			fprintf(gp, "plot [0:5][-1:1] \
					'%s' u 1:5 w l lw 2 lc 'red' title 'x', \
					'%s' u 1:6 w l lw 1.5 lc 'blue' dt (10,10) title 'y', \
					'%s' u 1:7 w l lw 1.5 lc 'green' dt (10,10,3,10) title 'z', \
					\n", ERROR_FILE, ERROR_FILE, ERROR_FILE);
		}
	}
	else
	{
		if (POS_F == 1)
		{
			if (LAG_F == 0)
				fprintf(gp, "plot [0:10][0:12]\
						'%s' u 1:4 w l lw 1.5  notitle, \
						\n", NEWTON_FILE);
			else
				fprintf(gp, "plot [0:10][0:12]\
						'%s' u 1:4 w l lw 1.5  notitle, \
						\n", LAG_FILE);
		}
		else if (POS_F == 0)
		{
			if (LAG_F == 0)
				fprintf(gp, "plot [0:5][-1:1]\
						'%s' u 1:5 w l lw 2 lc 'red' title 'φ', \
						'%s' u 1:6 w l lw 1.5 lc 'blue' dt (10,10) title 'θ', \
						'%s' u 1:7 w l lw 1.5 lc 'green' dt (10,10,3,10) title 'ψ', \
						\n", NEWTON_FILE, NEWTON_FILE, NEWTON_FILE);
			else
				fprintf(gp, "plot [0:5][-1:1] \
						'%s' u 1:5 w l lw 2 lc 'red' title 'φ', \
						'%s' u 1:6 w l lw 1.5 lc 'blue' dt (10,10) title 'θ', \
						'%s' u 1:7 w l lw 1.5 lc 'green' dt (10,10,3,10) title 'ψ', \
						\n", LAG_FILE, LAG_FILE, LAG_FILE);
		}
	}

	fflush(gp);
	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
