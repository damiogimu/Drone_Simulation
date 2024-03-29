#include "../sim.h"

size_t	search_newline(char *buf, size_t *i)
{
	*i = 0;
	while (buf[*i] != '\0')
	{
		if (buf[*i] == '\n')
			return (1);
		(*i)++;
	}
	return (0);
}

int		get_line_from_tmp(char **line, char **tmp)
{
	size_t	i;
	int		newline_flag;

	i = 0;
	newline_flag = search_newline(*tmp, &i);
	if (gnl_strjoin(line, *tmp, i) == -1)
		return (free_all(line, NULL, tmp));
	newline_flag == 1 ? ++i : i;
	if ((*tmp = (gnl_strdup(tmp, i))) == NULL)
		return (free_all(line, NULL, tmp));
	if ((*tmp)[0] == '\0')
		free_all(NULL, NULL, tmp);
	return (newline_flag);
}

ssize_t	gnl_read(int fd, char **line, char **tmp)
{
	size_t	i;
	ssize_t	rdsize;
	char	*buf;
	int		newline_flag;

	newline_flag = 0;
	while (newline_flag != 1)
	{
		if (!(buf = malloc(sizeof(char) * ((size_t)BUFFER_SIZE + 1))))
			return (free_all(line, NULL, NULL));
		if ((rdsize = read(fd, buf, (size_t)BUFFER_SIZE)) == -1)
			return (free_all(line, &buf, NULL));
		buf[rdsize] = '\0';
		if (rdsize == 0)
			return (read_final_line(&buf));
		newline_flag = search_newline(buf, &i);
		if ((gnl_strjoin(line, buf, i)) == -1)
			return (free_all(line, &buf, NULL));
		newline_flag == 1 ? ++i : i;
		if ((*tmp = gnl_strdup(&buf, i)) == NULL)
			return (free_all(line, &buf, NULL));
		if ((*tmp)[0] == '\0')
			free_all(NULL, NULL, tmp);
	}
	return (1);
}

int		get_next_line(int fd, char **line)
{
	static char	*tmp[256];
	ssize_t		result;
	int			is_finish;

	if (fd < 0 || 255 < fd || line == NULL || BUFFER_SIZE <= 0)
		return (-1);
	if (!(*line = malloc(sizeof(char) * 1)))
		return (-1);
	(*line)[0] = '\0';
	if (tmp[fd])
	{
		if ((is_finish = get_line_from_tmp(line, &tmp[fd])) == -1)
			return (-1);
		if (is_finish == 1)
			return (1);
	}
	if (((result = gnl_read(fd, line, &tmp[fd])) == -1))
		return (-1);
	return (result);
}
