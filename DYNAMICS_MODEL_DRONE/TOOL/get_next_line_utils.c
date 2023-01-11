/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   get_next_line_utils.c                              :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: dmogi <dmogi@student.42tokyo.jp>           +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/12/24 21:04:13 by dmogi             #+#    #+#             */
/*   Updated: 2022/11/24 17:30:00 by dmogi            ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../my_header.h"

size_t	gnl_strlen(char *str)
{
	size_t i;

	i = 0;
	while (str[i] != '\0')
		i++;
	return (i);
}

int		gnl_strjoin(char **line, char *src, size_t n)
{
	size_t	i;
	size_t	j;
	size_t	lsize;
	char	*save_line;

	i = 0;
	j = 0;
	lsize = gnl_strlen(*line);
	if ((save_line = gnl_strdup(line, 0)) == NULL)
		return (-1);
	if (!(*line = malloc(sizeof(char) * (lsize + n + 1))))
		return (free_all(&save_line, NULL, NULL));
	while (i < lsize)
	{
		(*line)[i] = save_line[i];
		i++;
	}
	free(save_line);
	while (j < n)
	{
		(*line)[i + j] = src[j];
		j++;
	}
	(*line)[i + j] = '\0';
	return (1);
}

char	*gnl_strdup(char **src, size_t nextl_p)
{
	size_t	i;
	size_t	new_size;
	char	*new;

	i = 0;
	new_size = gnl_strlen(&(*src)[nextl_p]);
	if (!(new = malloc(sizeof(char) * (new_size + 1))))
		return (NULL);
	if ((*src)[nextl_p] != '\0')
	{
		while (i < new_size)
		{
			new[i] = (*src)[nextl_p + i];
			i++;
		}
	}
	new[i] = '\0';
	free(*src);
	return (new);
}

int		free_all(char **line, char **buf, char **tmp)
{
	if (line != NULL)
	{
		free(*line);
		*line = NULL;
	}
	if (buf != NULL)
	{
		free(*buf);
		*buf = NULL;
	}
	if (tmp != NULL)
	{
		free(*tmp);
		*tmp = NULL;
	}
	return (-1);
}

int		read_final_line(char **buf)
{
	free(*buf);
	*buf = NULL;
	return (0);
}
