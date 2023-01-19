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

int setup(double ***state, double ***k, double **x)
{
	int i;

	*x = malloc(sizeof(double) * X_NUM);
	*state = malloc(sizeof(double *) * RK4_SIZE);
	*k = malloc(sizeof(double *) * RK4_SIZE);
	if (*state == NULL || *k == NULL || *x == NULL)
	{
		my_free(0, 0, *state);
		my_free(0, 0, *k);
		free(*x);
		return (-1);
	}
	i = 0;
	while(i < RK4_SIZE)
	{
		(*state)[i] = malloc(sizeof(double) * X_NUM);
		(*k)[i] = malloc(sizeof(double) * X_NUM);
		if ((*state)[i] == NULL || (*k)[i] == NULL)
		{
			my_free(1, i, *state);
			my_free(1, i, *k);
			free(*x);
			return (-1);
		}
		i++;
	}
	return (0);
}

void init(double **state, double **k, double *x, t_integral *intg, t_desire *des)
{
	int i, j;
	i = 0;
	while(i < RK4_SIZE)
	{
		des->ad_old[i] = des->dad_old[i] = 0.0;
		des->bd_old[i] = des->dbd_old[i] = 0.0;
		j = 0;
		while (j < INTG_N)
		{
			intg->e[i][j] = 0.0;
			intg->e_old[i][j] = 0.0;
			intg->area[i][j] = 0.0;
			j++;
		}
		j = 0;
		while (j < X_NUM)
		{
			state[i][j] = 0.0;
			k[i][j] = 0.0;
			j++;
		}
		i++;
	}
	x[0] = INIT_X;
	x[1] = INIT_Y;
	x[2] = INIT_Z;
	i = 3;
	while (i < X_NUM)
	{
		x[i] = 0.0;
		i++;
	}
}
