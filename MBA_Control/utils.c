#include "sim.h"

void my_free(int rev_f, int size, double **ptr)
{
	int i;

	if (rev_f == 1)
	{
		while (size >= 0)
		{
			free(ptr[size]);
			ptr[size] = NULL;
			size--;
		}
	}
	else
	{
		i = 0;
		while (i < size)
		{
			free(ptr[i]);
			ptr[i] = NULL;
			i++;
		}
	}
	free(ptr);
	ptr = NULL;
}

void my_fclose(int rev_f, int size, FILE **fd)
{
	int i;

	if (rev_f == 1)
	{
		while (size >= 0)
		{
			fclose(fd[size]);
			fd[size] = NULL;
			size--;
		}
	}
	else
	{
		i = 0;
		while (i < size)
		{
			fclose(fd[i]);
			fd[i] = NULL;
			i++;
		}
	}
}

int setup(double ***state, double ***k, FILE **fd, t_desire **des, t_rotor **rot, t_gif **gif)
{
	char output_files[FD_NUM][100]={PATH_FILE, ERROR_FILE, ANI_PATH_FILE, ANI_CABLE_FILE, ANI_XROTOR_FILE, ANI_YROTOR_FILE, DESIRE_FILE, ACCTIME_FILE, THRUST_FILE};
	int i;
	i = 0;
	while (i < FD_NUM)
	{
		if ((fd[i] = fopen(output_files[i], "w")) == NULL)
		{
			my_fclose(1, i, fd);
			return (-1);
		}
		i++;
	}
	*state = malloc(sizeof(double *) * RK4_SIZE);
	*k = malloc(sizeof(double *) * RK4_SIZE);
	*des = malloc(sizeof(t_desire));
	*rot = malloc(sizeof(t_rotor));
	*gif = malloc(sizeof(t_gif));
	if (*state == NULL || *k == NULL || *des == NULL || *rot == NULL || *gif == NULL)
	{
		my_free(0, 0, *state);
		my_free(0, 0, *k);
		free(*des);
		free(*rot);
		free(*gif);
		my_fclose(0, FD_NUM, fd);
		return (-1);
	}
	i = 0;
	while (i < RK4_SIZE)
	{
		(*state)[i] = malloc(sizeof(double) * X_NUM);
		(*k)[i] = malloc(sizeof(double) * X_NUM);
		if ((*state)[i] == NULL || (*k)[i] == NULL)
		{
			my_free(1, i, *state);
			my_free(1, i, *k);
			free(*des);
			free(*rot);
			free(*gif);
			my_fclose(0, FD_NUM, fd);
			return (-1);
		}
		i++;
	}
	return (0);
}

void init(double **k, t_desire *des, t_rotor *rot, t_gif *gif)
{
	int i, j;
	
	gif->accuracy = 1;
	while ((int)(dt*gif->accuracy) != 1)
		gif->accuracy *= 10;
	rot->coefi = (Cq/Ct)*(R/M_PI);
	des->ACC_T = 0.0;
	des->acc_t_f = 1;
	i = 0;
	while (i < RK4_SIZE)
	{
		des->ad_old[i] = 0.0;
		des->bd_old[i] = 0.0;
		des->dad_old[i] = 0.0;
		des->dbd_old[i] = 0.0;
		j = 0;
		while (j < X_NUM)
		{
			k[i][j] = 0.0;
			j++;
		}
		i++;
	}
}

void output_data1(double t, double *x, FILE **fd, t_gif *gif)
{
	double x_p, y_p, z_p;

	x_p = x[0] + r*cos(x[7])*sin(x[6]);
	y_p = x[1] + r*sin(x[7]);
	z_p = x[2] - r*cos(x[7])*cos(x[6]);

	fprintf(fd[0], "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ", \
					t,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12], \
					x[13],x[14],x[15],x_p,y_p,z_p);

	if ((int)(t*gif->accuracy) % (int)(GIF_ITV*gif->accuracy) == 0)  // GIT_ITV 毎に出力
	{
		int f = 1;
		if (f == 1)
		{
			fprintf(fd[3], "%f %f %f %f\n", t, x[0], x[1], x[2]);
			fprintf(fd[3], "%f %f %f %f\n", t, x_p, y_p, z_p);
			fprintf(fd[3], "\n\n");
			gif->l_x1 = x[0]+(SCALE*l)*cos(x[3])*cos(x[5]);
			gif->l_y1 = x[1]+(SCALE*l)*sin(x[4]);
			gif->l_z1 = x[2]-(SCALE*l)*tan(x[4]);
			fprintf(fd[4], "%f %f %f %f\n", t, gif->l_x1, gif->l_y1, gif->l_z1);
			gif->l_x2 = x[0]-(SCALE*l)*cos(x[5])*cos(x[5]);
			gif->l_y2 = x[1]-(SCALE*l*cos(x[5])*sin(x[5]));
			gif->l_z2 = x[2] + (SCALE*l)*tan(x[4]);
			fprintf(fd[4], "%f %f %f %f\n", t, gif->l_x2, gif->l_y2, gif->l_z2);
			fprintf(fd[4], "\n\n");
			gif->l_x1 = x[0]+(SCALE*l*cos(x[5])*sin(x[5]));
			gif->l_y1 = x[1]-(SCALE*l)*cos(x[5])*cos(x[5]);
			gif->l_z1 = x[2] - (SCALE*l)*tan(x[3]);
			fprintf(fd[5], "%f %f %f %f\n", t, gif->l_x1, gif->l_y1, gif->l_z1);
			gif->l_x2 = x[0]-(SCALE*l*sin(x[5])*cos(x[5]));
			gif->l_y2 = x[1]+(SCALE*l)*cos(x[5])*cos(x[5]);
			gif->l_z2 = x[2] + (SCALE*l)*tan(x[3]);
			fprintf(fd[5], "%f %f %f %f\n", t, gif->l_x2, gif->l_y2, gif->l_z2);
			fprintf(fd[5], "\n\n");
			fprintf(fd[2], "%f %f %f %f %f %f %f %f %f %f", t,x[0],x[1],x[2],x[3],x[4],x[5],x_p,y_p,z_p);
		}
		else
		{
			fprintf(fd[3], "%f %f %f %f\n", t, x[0], x[1], x[2]);
			fprintf(fd[3], "%f %f %f %f\n", t, x_p, y_p, z_p);
			fprintf(fd[3], "\n\n");
			gif->l_x1 = x[0]+(SCALE*l)*cos(x[5])*cos(x[5]);
			gif->l_y1 = x[1]+(SCALE*l)*sin(x[5])*cos(x[5]);
			gif->l_z1 = x[2]-(SCALE*l)*tan(x[4]);
			fprintf(fd[4], "%f %f %f %f\n", t, gif->l_x1, gif->l_y1, gif->l_z1);
			gif->l_x2 = x[0]-(SCALE*l)*cos(x[5])*cos(x[5]);
			gif->l_y2 = x[1]-(SCALE*l*cos(x[5])*sin(x[5]));
			gif->l_z2 = x[2] + (SCALE*l)*tan(x[4]);
			fprintf(fd[4], "%f %f %f %f\n", t, gif->l_x2, gif->l_y2, gif->l_z2);
			fprintf(fd[4], "\n\n");
			gif->l_x1 = x[0]+(SCALE*l*cos(x[5])*sin(x[5]));
			gif->l_y1 = x[1]-(SCALE*l)*cos(x[5])*cos(x[5]);
			gif->l_z1 = x[2] - (SCALE*l)*tan(x[3]);
			fprintf(fd[5], "%f %f %f %f\n", t, gif->l_x1, gif->l_y1, gif->l_z1);
			gif->l_x2 = x[0]-(SCALE*l*sin(x[5])*cos(x[5]));
			gif->l_y2 = x[1]+(SCALE*l)*cos(x[5])*cos(x[5]);
			gif->l_z2 = x[2] + (SCALE*l)*tan(x[3]);
			fprintf(fd[5], "%f %f %f %f\n", t, gif->l_x2, gif->l_y2, gif->l_z2);
			fprintf(fd[5], "\n\n");
			fprintf(fd[2], "%f %f %f %f %f %f %f %f %f %f", t,x[0],x[1],x[2],x[3],x[4],x[5],x_p,y_p,z_p);
		}
	}
}

void output_data2(double t, double *x, double **k, FILE **fd, t_desire *des, t_rotor *rot, t_gif *gif)	
{
	double ddx, ddy, ddz;
		
	ddx = (k[0][8] + 2.0*k[1][8] + 2.0*k[2][8] + k[3][8])/6.0;
	ddy = (k[0][9] + 2.0*k[1][9] + 2.0*k[2][9] + k[3][9])/6.0;
	ddz = (k[0][10] + 2.0*k[1][10] + 2.0*k[2][10] + k[3][10])/6.0;
	fprintf(fd[0], "%f %f %f\n", ddx, ddy, ddz);
	fprintf(fd[6], "%f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, des->ddxd, des->ddyd, des->ddzd, des->dxd, des->dyd, des->dzd, des->xd, des->yd, des->zd, des->ad, des->bd, des->cd);
	fprintf(fd[8], "%f %f %f %f %f %f\n", t, rot->T1, rot->T2, rot->T3, rot->T4, rot->T1+rot->T2+rot->T3+rot->T4);
	if ((int)(t*gif->accuracy) % (int)(GIF_ITV*gif->accuracy) == 0)
		fprintf(fd[2], " %f %f %f\n", des->xd, des->yd, des->zd);
	fprintf(fd[1], "%f %f %f %f ", t, (des->xd-x[0]), (des->yd-x[1]), (des->zd-x[2]));
	if (des->acc_t_f == 1)
	{
		fprintf(fd[7], "%f\n", des->ACC_T);
		des->acc_t_f = 0;
	}
}
