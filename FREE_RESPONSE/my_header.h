#ifndef MY_HEADER_H
# define MY_HEADER_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define amax 1.0	// 最大加速度
#define vmax 0.5	// 終端速度
#define Xt 2.0		// 終端位置 x
#define Yt 0.0		// 終端位置 y

#define TIME 7.0

#define INIT_X 0.0
#define INIT_Y 0.0
#define INIT_Z 2.0

#define RK4_SIZE 4
#define X_NUM 16
#define dt 1.0e-3

typedef struct	s_rotor
{
	double	T1, T2, T3, T4;
	double	coefi;
	double	l_x1, l_x2;
	double	l_y1, l_y2;
	double	l_z1, l_z2;
}				t_rotor;

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
#define FD_NUM 5
#define PATH_FILE "DATA/result_line_path"
#define ANI_PATH_FILE "DATA/path_DATA"
#define ANI_CABLE_FILE "DATA/cable_DATA"
#define ANI_XROTOR_FILE "DATA/xrotor_DATA"
#define ANI_YROTOR_FILE "DATA/yrotor_DATA"
// ------------------------ //

int		setup(double ***state, double ***k, FILE **fd);
void	init(double **k, t_rotor *rot);
void	output_data1(double t, double *x, FILE **fd, t_rotor *rot);
void	my_free(int rev_f, int size, double **ptr);
void	my_fclose(int rev_f, int size, FILE **fd);

int		func(double t, double *x, double *k, int i, t_rotor *rot);

#endif
