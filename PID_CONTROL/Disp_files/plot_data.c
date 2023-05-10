#include "../sim.h"

#define X_MIN 0.0
#define X_ITV 1.0

#define PDF_NUM 4
#define ACC_F "Graphs/acc.pdf"
#define VELO_F "Graphs/velo.pdf"
#define POS_F "Graphs/pos.pdf"
#define ATTI_F "Graphs/attitude.pdf"
#define LOAD_F "Graphs/load_gamma.pdf"

#define LABEL_ACC "Acceleration [m/s^2]"
#define LABEL_VELO "Velocity [m/s]"
#define LABEL_POS "Postion [m]"
#define LABEL_ANGLE "Angle [rad]"

#define ACC_T_FILE "Data/acc_time"

#define COL_X_ACC 21
#define COL_X_VELO 10
#define COL_X_POS 2
#define COL_THETA 6

#define COL_X_ACC_DES 2
#define COL_X_VELO_DES 5
#define COL_X_POS_DES 8
#define COL_THETA_DES 12

int main(void)
{
	FILE *gp;
	char *line;
	char pdf_file_name[PDF_NUM][50]={ACC_F,VELO_F,POS_F,ATTI_F};
	char label_name[PDF_NUM][30]={LABEL_ACC,LABEL_VELO,LABEL_POS,LABEL_ANGLE};
	char data_title[PDF_NUM][10]={"x acc", "x velo", "x pos", "theta"};
	int i;
	int fd;
	int ref_data_col[PDF_NUM]={COL_X_ACC,COL_X_VELO,COL_X_POS,COL_THETA};
	int des_data_col[PDF_NUM]={COL_X_ACC_DES,COL_X_VELO_DES,COL_X_POS_DES,COL_THETA_DES}; 
	double acc_t;
	double y_max[PDF_NUM+1]={20, vmax+1.0, Xt+1.0, 0.3, 0.6};
	double y_min[PDF_NUM+1]={-(20), -(vmax+1.0), -(Xt+1.0), -0.3, -0.6};

	if ((gp = popen("/usr/local/bin/gnuplot", "w")) == NULL || (fd = open("DATA/acc_time", O_RDONLY)) == -1)
	{
		fprintf(stderr, "OPEN ERROR:file does not exist\n");
		return (-1);
	}
	get_next_line(fd, &line);
	acc_t = atof(line) + Z_RISE_T;


	free(line);
	close(fd);

	fprintf(gp, "set size ratio 0.6\n");
	fprintf(gp, "set grid\n");
	fprintf(gp, "set xtics %f,%f,%f\n", X_MIN, X_ITV, TIME);
	fprintf(gp, "set xlabel 'Time [s]'\n");
	fprintf(gp, "set border lw 3\n");

	fprintf(gp, "set term pdf size 6in, 4in font 'Times New Roman, 26'\n");
	fprintf(gp, "set key font 'Times New Roman, 22'\n");


	i = 0;
	while (i < PDF_NUM)
	{
		if (i == 0)
			fprintf(gp, "set arrow 1 from %f,%f to %f,%f nohead dt (5, 15)\n", acc_t, y_min[i], acc_t, y_max[i]);
		else
			fprintf(gp, "unset arrow\n");
		fprintf(gp, "set ylabel '%s'\n", label_name[i]);
		fprintf(gp, "set output '%s'\n", pdf_file_name[i]);
		fprintf(gp, "plot [%f:%f][%f:%f]\
				'Data/result' u 1:%d w l lc 'blue' title '%s', \
				'Data/desire' u 1:%d w l lc 'red' dt(10,10) title 'desire', \
				\n", X_MIN, TIME, y_min[i], y_max[i], ref_data_col[i], data_title[i], des_data_col[i]);
		i++;
	}
	fprintf(gp, "set arrow 1 from %f,%f to %f,%f nohead dt (5, 15)\n", acc_t, y_min[i], acc_t, y_max[i]);
	fprintf(gp, "set ylabel '%s'\n", LABEL_ANGLE);
	fprintf(gp, "set output '%s'\n", LOAD_F);
	fprintf(gp, "plot [%f:%f][%f:%f]'Data/result' u 1:8 w l lc 'blue' title 'gamma'\n", X_MIN, TIME, y_min[i], y_max[i]);

	fflush(gp);		// gpとのストリーム上のバッファのデータを放出する
	fprintf(gp, "exit\n");
	pclose(gp);
	return (0);
}
