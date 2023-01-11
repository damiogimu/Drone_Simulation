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
