#include "my_header.h"

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

int setup(double ***state, double ***k, FILE **fd)
{
	char output_files[FD_NUM][100]={PATH_FILE, ERROR_FILE, ANI_PATH_FILE, ANI_CABLE_FILE, ANI_XROTOR_FILE, ANI_YROTOR_FILE, DESIRE_FILE};
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
	if (*state == NULL || *k == NULL)
	{
		my_free(0, 0, *state);
		my_free(0, 0, *k);
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
			return (-1);
		}
		i++;
	}
	return (0);
}

void init(double **k, t_rotor *rot, t_desire *des)
{
	int i, j;

	rot->coefi = (Cq/Ct)*(R/M_PI);
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

void output_data1(double t, double *x, FILE **fd, t_rotor *rot)
{
	double x_p, y_p, z_p;

	x_p = x[0] + r*cos(x[7])*sin(x[6]);
	y_p = x[1] + r*sin(x[7]);
	z_p = x[2] - r*cos(x[7])*cos(x[6]);
	fprintf(fd[0], "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", \
					t,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12], \
					x[13],x[14],x[15],x_p,y_p,z_p);

	if ((int)((t*1000000.0)+0.5)%50000 == 0)
	{
		fprintf(fd[3], "%f %f %f %f\n", t, x[0], x[1], x[2]);
		fprintf(fd[3], "%f %f %f %f\n", t, x_p, y_p, z_p);
		fprintf(fd[3], "\n\n");
		rot->l_x1 = x[0]+(SCALE*l)*cos(x[5])*cos(x[5]);
		rot->l_y1 = x[1]+(SCALE*l)*sin(x[5])*cos(x[5]);
		rot->l_z1 = x[2] -(SCALE*l)*tan(x[4]);
		fprintf(fd[4], "%f %f %f %f\n", t, rot->l_x1, rot->l_y1, rot->l_z1);
		rot->l_x2 = x[0]-(SCALE*l)*cos(x[5])*cos(x[5]);
		rot->l_y2 = x[1]-(SCALE*l*cos(x[5])*sin(x[5]));
		rot->l_z2 = x[2] + (SCALE*l)*tan(x[4]);
		fprintf(fd[4], "%f %f %f %f\n", t, rot->l_x2, rot->l_y2, rot->l_z2);
		fprintf(fd[4], "\n\n");
		rot->l_x1 = x[0]+(SCALE*l*cos(x[5])*sin(x[5]));
		rot->l_y1 = x[1]-(SCALE*l)*cos(x[5])*cos(x[5]);
		rot->l_z1 = x[2] - (SCALE*l)*tan(x[3]);
		fprintf(fd[5], "%f %f %f %f\n", t, rot->l_x1, rot->l_y1, rot->l_z1);
		rot->l_x2 = x[0]-(SCALE*l*sin(x[5])*cos(x[5]));
		rot->l_y2 = x[1]+(SCALE*l)*cos(x[5])*cos(x[5]);
		rot->l_z2 = x[2] + (SCALE*l)*tan(x[3]);
		fprintf(fd[5], "%f %f %f %f\n", t, rot->l_x2, rot->l_y2, rot->l_z2);
		fprintf(fd[5], "\n\n");
		fprintf(fd[2], "%f %f %f %f %f %f %f %f %f %f", t,x[0],x[1],x[2],x[3],x[4],x[5],x_p,y_p,z_p);
	}
}

void output_data2(double t, double *x, double **k, FILE **fd, t_rotor *rot, t_desire *des)	
{
	double ddx, ddy, ddz;
		
	ddx = (k[0][8] + 2.0*k[1][8] + 2.0*k[2][8] + k[3][8])/6.0;
	ddy = (k[0][9] + 2.0*k[1][9] + 2.0*k[2][9] + k[3][9])/6.0;
	ddz = (k[0][10] + 2.0*k[1][10] + 2.0*k[2][10] + k[3][10])/6.0;
	fprintf(fd[6], "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, des->ddxd, des->ddyd, des->ddzd, ddx, ddy, ddz, des->dxd, des->dyd, des->dzd, des->xd, des->yd, des->zd, des->ad, des->bd, des->cd);
	if ((int)((t*1000000.0)+0.5)%50000==0)
		fprintf(fd[2], " %f %f %f\n", des->xd, des->yd, des->zd);
	fprintf(fd[1], "%f %f %f %f ", t, (des->xd-x[0]), (des->yd-x[1]), (des->zd-x[2]));
}
