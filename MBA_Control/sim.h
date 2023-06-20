#ifndef MY_HEADER_H
# define MY_HEADER_H

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define amax 1.0
#define vmax 0.5
#define Xt 2.0
#define Yt 0.0

#define TIME 14.0
#define Z_RISE_T 1.0
#define dt 1.0e-3

#define INIT_X 0.0
#define INIT_Y 0.0
#define INIT_Z 2.0

#define RK4_SIZE 4
#define X_NUM 16
#define GIF_ITV 0.05

typedef struct	s_desire
{
	int		acc_t_f;
	double	ACC_T;
	double	xd, dxd, ddxd;
	double	yd, dyd, ddyd;
	double	zd, dzd, ddzd;
	double	ad, dad, ddad, ad_old[RK4_SIZE], dad_old[RK4_SIZE];
	double	bd, dbd, ddbd, bd_old[RK4_SIZE], dbd_old[RK4_SIZE];
	double	cd, dcd, ddcd;
}				t_desire;

typedef struct	s_rotor
{
	double T1, T2, T3, T4;
	double T1_o, T2_o, T3_o, T4_o;
	double coefi;
}				t_rotor;

typedef struct	s_gif
{
	int	accuracy;
	double l_x1, l_x2;
	double l_y1, l_y2;
	double l_z1, l_z2;
}				t_gif;

// --- CONTRLER GAIN VALS --- //
#define Kpt 70
#define Kdt 50
#define Kpr 70
#define Kdr 70
// ------------------------- //

// ---- PARAMETERS ---- //
#define m 0.65
#define mp 0.3
#define Ixx 7.5e-3
#define Iyy 7.5e-3
#define Izz 1.3e-2

#define SCALE 1.0
#define l 0.232			// プロペラ中心と質量中心の距離
#define r 1.0 			// ケーブルの長さ
#define Ct 0.07428
#define Cq 0.10724
#define Kt 10.0e-15
#define Kr 10.0e-15

#define g 9.81
#define p 1.29			// 空気の密度 [kg/m^3]
#define R 0.18			// プロペラの半径 [m]
// ------------------ //

// --- OUTPUT DATA FILES --- //
#define FD_NUM 9
#define PATH_FILE "DATA/result"
#define ERROR_FILE "DATA/error"
#define ANI_PATH_FILE "DATA/GIF_path"
#define ANI_CABLE_FILE "DATA/GIF_cable"
#define ANI_XROTOR_FILE "DATA/GIF_xrotor"
#define ANI_YROTOR_FILE "DATA/GIF_yrotor"
#define DESIRE_FILE "DATA/desire"
#define ACCTIME_FILE "DATA/acc_time"
#define THRUST_FILE "DATA/rotor_thrust"
// ------------------------ //

void	any_traj(double t, t_desire *des);
void	noncontrol_traj_vt(double t, t_desire *des);
void	noncontrol_traj_xt(double t, t_desire *des);

int		setup(double ***state, double ***k, FILE **fd, t_desire **des, t_rotor **rot, t_gif **gif);
void	init(double **k, t_desire *des, t_rotor *rot, t_gif *gif);
void	output_data1(double t, double *x, FILE **fd, t_gif *gif);
void	output_data2(double t, double *x, double **k, FILE **fd, t_desire *des, t_rotor *rot, t_gif *gif);

void	my_free(int rev_f, int size, double **ptr);
void	my_fclose(int rev_f, int size, FILE **fd);
int		func(double t, double *x, double *k, int i, t_rotor *rot, t_desire *des);
int		get_next_line(int fd, char **line);

#ifndef BUFFER_SIZE
# define BUFFER_SIZE 100
#endif
size_t	gnl_strlen(char *str);
int		gnl_strjoin(char **line, char *src, size_t n);
char	*gnl_strdup(char **src, size_t newl_p);
int		free_all(char **line, char **buf, char **tmp);
size_t	search_newline(char *buf, size_t *i);
int		get_line_from_tmp(char **line, char **tmp);
int		read_final_line(char **buf);
ssize_t	gnl_read(int fd, char **line, char **tmp);
#endif
