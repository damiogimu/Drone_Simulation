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
	char output_files[FD_NUM][100]={PATH_FILE, ANI_PATH_FILE, ANI_CABLE_FILE, ANI_XROTOR_FILE, ANI_YROTOR_FILE};
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
			my_fclose(0, FD_NUM, fd);
			return (-1);
		}
		i++;
	}
	return (0);
}

void init(double **k, t_rotor *rot)
{
	int i, j;

	rot->coefi = (Cq/Ct)*(R/M_PI);
	i = 0;
	while (i < RK4_SIZE)
	{
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
		fprintf(fd[2], "%f %f %f %f\n", t, x[0], x[1], x[2]);
		fprintf(fd[2], "%f %f %f %f\n", t, x_p, y_p, z_p);
		fprintf(fd[2], "\n\n");
		rot->l_x1 = x[0]+(SCALE*l)*cos(x[5])*cos(x[5]);
		rot->l_y1 = x[1]+(SCALE*l)*sin(x[5])*cos(x[5]);
		rot->l_z1 = x[2] -(SCALE*l)*tan(x[4]);
		fprintf(fd[3], "%f %f %f %f\n", t, rot->l_x1, rot->l_y1, rot->l_z1);
		rot->l_x2 = x[0]-(SCALE*l)*cos(x[5])*cos(x[5]);
		rot->l_y2 = x[1]-(SCALE*l*cos(x[5])*sin(x[5]));
		rot->l_z2 = x[2] + (SCALE*l)*tan(x[4]);
		fprintf(fd[3], "%f %f %f %f\n", t, rot->l_x2, rot->l_y2, rot->l_z2);
		fprintf(fd[3], "\n\n");
		rot->l_x1 = x[0]+(SCALE*l*cos(x[5])*sin(x[5]));
		rot->l_y1 = x[1]-(SCALE*l)*cos(x[5])*cos(x[5]);
		rot->l_z1 = x[2] - (SCALE*l)*tan(x[3]);
		fprintf(fd[4], "%f %f %f %f\n", t, rot->l_x1, rot->l_y1, rot->l_z1);
		rot->l_x2 = x[0]-(SCALE*l*sin(x[5])*cos(x[5]));
		rot->l_y2 = x[1]+(SCALE*l)*cos(x[5])*cos(x[5]);
		rot->l_z2 = x[2] + (SCALE*l)*tan(x[3]);
		fprintf(fd[4], "%f %f %f %f\n", t, rot->l_x2, rot->l_y2, rot->l_z2);
		fprintf(fd[4], "\n\n");
		fprintf(fd[1], "%f %f %f %f %f %f %f %f %f %f\n", t,x[0],x[1],x[2],x[3],x[4],x[5],x_p,y_p,z_p);
	}
}
