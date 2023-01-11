#include "my_header.h"

#ifndef POS_F
# define POS_F 0
#endif

#define NEWTON_FILE "./RESULT/result_newton"
#define LAG_FILE "./RESULT/result_lag"

#define X_MIN 0.0
#define X_MAX 5.0

#define Y_MIN 0.0
#define Y_MAX 1.0

#define STATE_NUM 6

int main(void)
{
	FILE *gp;
	char buf[5];
	char *line1,*line2; 
	int i;
	int fd1, fd2;
	int gnl_f1=1, gnl_f2=1;
	int start1=0, end1=0, start2=0, end2=0;
	double t = 0.0;
	double q_lag[STATE_NUM];
	double q_nt[STATE_NUM];

// --- fd setup --- //
	if ((gp = popen("/usr/local/bin/gnuplot", "w")) == NULL)
	{
		fprintf(stderr, "Gnuplot exe_file couldn't open\n");
		return (-1);
	}
	if ((fd1 = open(LAG_FILE, O_RDONLY)) == -1)
	{
		fprintf(stderr, "lagrange_result_file could't open\n'");
		pclose(gp);
		return (-1);
	}
	if ((fd2 = open(NEWTON_FILE, O_RDONLY)) == -1)
	{
		fprintf(stderr, "newton_result_file could't open\n'");
		pclose(gp);
		close(fd1);
		return (-1);
	}
// ---------------- //

// --- calc diff --- //

	while (1)
	{
		start1 = start2 = 0;
		end1 = end2 = 0;
		gnl_f1 = get_next_line(fd1, &line1);
		gnl_f2 = get_next_line(fd2, &line2);
		if (gnl_f1 == 0 || gnl_f2 == 0)
			break;
		start1 = skip_space(&line1[start1]);
		end1 = skip_space(&line1[start1]) + skip_num(&line1[start1]);
		t = atof(&line1[start1]);
		end2 = skip_space(&line2[start2]) + skip_num(&line2[start2]);
		i = 0;
		while (i < STATE_NUM)
		{
			start1 = end1 + skip_space(&line1[end1]);
			end1 = start1 + skip_num(&line1[start1]);
			start2 = end2 + skip_space(&line2[end2]);
			end2 = start2 + skip_num(&line2[start2]);
			q_lag[i] = atof(&line1[start1]);
			q_nt[i] = atof(&line2[start2]);
			line1[end1] = '\0';
			line2[end2] = '\0';
			end1++;
			end2++;
			i++;
		}
		if (gnl_f1 != 0 && gnl_f2 != 0)
			printf("%f %f %f %f %f %f %f\n", t, q_lag[0]-q_nt[0], q_lag[1]-q_nt[1], q_lag[2]-q_nt[2], q_lag[3]-q_nt[3], q_lag[4]-q_nt[4], q_lag[5]-q_nt[5]);
		free(line1);
		free(line2);
		line1 = line2 = NULL;
	}

// ----------------- //

/*
// --- gnuplot setup --- //
	fprintf(gp, "set terminal qt font 'Times New Roman, 20'\n");
	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set xtics \n");
	fprintf(gp, "set ytics \n");
	fprintf(gp, "set xlabel 'T [s]'\n");
	if (POS_F == 1)
	{
		fprintf(gp, "set ylabel 'Position [m]'\n");
		fprintf(gp, "set key at 7,4.5\n");
		fprintf(gp, "set xrange[0:10]\n");
		fprintf(gp, "set yrange[0:5]\n");
	}
	else
	{
		fprintf(gp, "set ylabel 'angle [rad]'\n");
		fprintf(gp, "set key at 4,2.8\n");
		fprintf(gp, "set xrange[0:5]\n");
		fprintf(gp, "set yrange[0:3]\n");
	}
	fprintf(gp, "set border lw 3\n");
	fprintf(gp, "set key font ',16'\n");
// --------------------- //

// --- plot diff lag_newton --- //
	if (POS_F == 1)
	{
		fprintf(gp, "plot \
				'%s' u 1:2 w l lw 2 lc 'red' title 'x erorr', \
				'%s' u 1:3 w l lw 1.5 lc 'blue' dt (10,10) title 'y error', \
				'%s' u 1:4 w l lw 1.5 lc 'dark-green' dt (10,10,3,10) title 'z error', \
				\n", POS_DIFF, POS_DIFF, POS_DIFF);
	}
	else
	{
		fprintf(gp, "plot \
				'%s' u 1:5 w l lw 2 lc 'red' title 'φ error', \
				'%s' u 1:6 w l lw 1.5 lc 'blue' dt (10,10) title 'θ error', \
				'%s' u 1:7 w l lw 1.5 lc 'dark-green' dt (10,10,3,10) title 'ψ error', \
				\n", ANGLE_DIFF, ANGLE_DIFF, ANGLE_DIFF);
	}
// ---------------------------- //
*/
	fflush(gp);
//	read(0, buf, 1);
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);	
}
