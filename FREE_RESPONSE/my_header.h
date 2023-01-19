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
# define TIME 10.0
#endif

#define INIT_X 0.0
#define INIT_Y 0.0
#define INIT_Z 1.0

#define RK4_SIZE 4
#define X_NUM 16
#define dt 5.0e-4

typedef struct	s_rotor
{
	double T1, T2, T3, T4;
	double coefi;
	double l_x1, l_x2;
	double l_y1, l_y2;
	double l_z1, l_z2;
	
}				t_rotor;

// ---- Parameter ---- //
#define m 2.0
#define mp 0.1
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

#define FD_NUM 5
#define DATA_FILE1 "DATA/result_lag"
#define DATA_FILE2 "../ANIMATION_DYNAMICS/DATA/cable_DATA"
#define DATA_FILE3 "../ANIMATION_DYNAMICS/DATA/path_DATA"
#define DATA_FILE4 "../ANIMATION_DYNAMICS/DATA/xrotor_DATA"
#define DATA_FILE5 "../ANIMATION_DYNAMICS/DATA/yrotor_DATA"

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
