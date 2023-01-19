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
	char output_files[FD_NUM][100]={PATH_FILE, ERROR_FILE, ANI_PATH_FILE, ANI_CABLE_FILE, ANI_XROTOR_FILE, ANI_YROTOR_FILE, THRUST_TORQU, DESIRE_FILE};
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
