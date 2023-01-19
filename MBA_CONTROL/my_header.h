#ifndef MY_HEADER_H
# define MY_HEADER_H

#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define amax 3.0
#define vmax 3.0
#define Xt 2.0
#define Yt 0.0

#define TIME 8.0
#define Z_RISE_T 0.0
#define INIT_X 0.0
#define INIT_Y 0.0
#define INIT_Z 2.0

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
	double T1, T2, T3, T4;
	double coefi;
	double l_x1, l_x2;
	double l_y1, l_y2;
	double l_z1, l_z2;
}				t_rotor;

// --- CONTROLER GAIN VALS --- //
#define Kpt 15
#define Kdt 15
#define Kpr 15
#define Kdr 55
// ------------------------- //

// ---- PARAMETERS ---- //
#define m 0.65
#define mp 0.3
#define Ixx 7.5e-3
#define Iyy 7.5e-3
#define Izz 1.3e-2

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
#define FD_NUM 8
#define PATH_FILE "DATA/result_line_path"
#define ERROR_FILE "DATA/result_line_error"
#define ANI_PATH_FILE "../ANIMATION_MBA/DATA/path_DATA"
#define ANI_CABLE_FILE "../ANIMATION_MBA/DATA/cable_DATA"
#define ANI_XROTOR_FILE "../ANIMATION_MBA/DATA/xrotor_DATA"
#define ANI_YROTOR_FILE "../ANIMATION_MBA/DATA/yrotor_DATA"
#define THRUST_TORQU "DATA/rotor_thrust_and_torqu"
#define DESIRE_FILE "DATA/desire_path"
// ------------------------ //

void	any_traj(double t, t_desire *des);
void	noncontrol_traj1(double t, t_desire *des);
void	noncontrol_traj2(double t, t_desire *des);
void	controled_traj1(double t, double f, t_desire *des);
void	controled_traj2(double t, double f, t_desire *des);

int		setup(double ***state, double ***k, FILE **fd);
void	init(double **k, t_rotor *rot, t_desire *des);
void	my_free(int rev_f, int size, double **ptr);
void	my_fclose(int rev_f, int size, FILE **fd);

int		func(double t, double *x, double *k, int i, t_rotor *rot, t_desire *des);

#endif
