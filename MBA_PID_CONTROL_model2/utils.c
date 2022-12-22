/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   utils.c                                            :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: dmogi <mogi.t218140@gmail.com>             +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2022/12/22 18:00:02 by dmogi             #+#    #+#             */
/*   Updated: 2022/12/22 18:23:28 by dmogi            ###   ########.jp       */
/*                                                                            */
/* ************************************************************************** */

#include "my_header.h"

void my_free(int size, double **ptr)
{
	int i;
	i = 0;
	while (i < size)
	{
		free(ptr[i]);
		ptr[i] = NULL;
		i++;
	}
	free(ptr);
	ptr = NULL;
}

void my_fclose(int size, FILE **fd)
{
	int i;
	i = 0;
	while (i < size)
	{
		fclose(fd[i]);
		fd[i] = NULL;
		i++;
	}	
}
