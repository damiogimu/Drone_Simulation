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

#define amax 1.0
#define vmax 0.5
#define Xt 2.0
#define Yt 2.0

#define TIME 5.0
#define INIT_X 0.0
#define INIT_Y 0.0
#define INIT_Z 2.0

#define Z_RISE_T 1.0
#define RK4_SIZE 4
#define X_NUM 16
#define dt 1.0e-3

typedef struct	s_desire
{
	double	xd, dxd, ddxd;
	double	yd, dyd, ddyd;
	double	zd, dzd, ddzd;
	double	ad, dad, ddad, ad_old[RK4_SIZE], dad_old[RK4_SIZE];
	double	bd, dbd, ddbd, bd_old[RK4_SIZE], dbd_old[RK4_SIZE];
	double	cd, dcd, ddcd;
}				t_desire;

typedef struct	s_rotor 
{
	double	T1, T2, T3, T4;
	double	coefi;
	double	l_x1, l_x2;
	double	l_y1, l_y2;
	double	l_z1, l_z2;
}				t_rotor;

typedef struct	s_gain
{
	int Kpt, Kdt;
	int Kpr, Kdr;
}				t_gain;

// --- GAIN_RANGE --- //
#define MIN_Kpt_GAIN 40
#define MAX_Kpt_GAIN 80
#define Kpt_ITV 10
#define MIN_Kdt_GAIN 10
#define MAX_Kdt_GAIN 80
#define Kdt_ITV 10
#define MIN_Kpr_GAIN 60
#define MAX_Kpr_GAIN 80
#define Kpr_ITV 10
#define MIN_Kdr_GAIN 70
#define MAX_Kdr_GAIN 80
#define Kdr_ITV 10
// ------------------- //

// ---- Parameter ---- //
#define m 0.65
#define mp 0.3
#define Ixx 7.5e-3
#define Iyy 7.5e-3
#define Izz 1.3e-2

#define l 0.232			// プロペラ中心と質量中心の距離
#define r 1.0
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

void	my_free(int rev_f, int size, double **ptr);
void	my_fclose(int rev_f, int size, FILE **fd);
void	set_desire_path1(double t, double f, t_desire *des);
void	set_desire_path2(double t, double f, t_desire *des);
void	init(double **state, double **k, double *x, t_desire *des);
int		setup(double ***state, double ***k, double **x);
int		func(double t, double *x, double *k, int i, t_desire *des, t_gain *gain);

#endif
