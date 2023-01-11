#ifndef MY_HEADER_H
# define MY_HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef TIME
# define TIME 6.0
#endif

#ifndef X_POS
# define X_POS 0.0
#endif
#ifndef Y_POS
# define Y_POS 0.0
#endif
#ifndef Z_POS
# define Z_POS 2.0
#endif

#define PHI 0.0 //
#define THETA 0.0 // 加速度 0.01 の時の角度
#define PSI 0.0

#define RK4_SIZE 4
#define X_NUM 16
#define dt 5.0e-4
#define PI 3.1415

// ---- Parameter ---- //
#define m 0.65
#define mp 0.3 // 0.3
#define Ixx 7.5e-3
#define Iyy 7.5e-3
#define Izz 1.3e-2

#define l 0.232			// プロペラ中心と質量中心の距離
#define r 1.0 // 1.0 			
#define Ct 0.07428		
#define Cq 0.10724		
#define Kt 10.0e-15		
#define Kr 10.0e-15		

#define g 9.81			
#define p 1.29			// 空気の密度 [kg/m^3]
#define R 0.18			// プロペラの半径 [m]
// ----------------- //

// ---- get_next_line ---- //
#define BUFFER_SIZE 100
size_t	gnl_strlen(char *str);
int		gnl_strjoin(char **line, char *src, size_t n);
char	*gnl_strdup(char **src, size_t newl_p);
int		free_all(char **line, char **buf, char **tmp);
size_t	search_newline(char *buf, size_t *i);
int		get_line_from_tmp(char **line, char **tmp);
int		read_final_line(char **buf);
ssize_t	gnl_read(int fd, char **line, char **tmp);
int		get_next_line(int fd, char **line);
int 	skip_num(char *line);
int		skip_space(char *line);
// ----------------------- //

#endif
