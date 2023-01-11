#ifndef MY_HEADER_H
# define MY_HEADER_H

#include<unistd.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define INIT_X 0.0
#define INIT_Y 0.0
#define INIT_Z 1.0
#define TIME 30.0

#define RK4_SIZE 4
#define PID_GAIN_NUM 6
#define X_NUM 16
#define dt 1.0e-3

typedef struct	s_integral
{
	double e[RK4_SIZE][PID_GAIN_NUM];
	double e_old[RK4_SIZE][PID_GAIN_NUM];
	double area[RK4_SIZE][PID_GAIN_NUM];
}				t_integral;

typedef struct	s_desire
{
	double xd, dxd, ddxd;
	double yd, dyd, ddyd;
	double zd, dzd, ddzd;
	double ad, dad, ddad;
	double bd, dbd, ddbd;
	double cd, dcd, ddcd;
}				t_desire;

typedef struct	s_rotor
{
	double T1;
	double T2;
	double T3;
	double T4;
	double coefi;
	double l_x1, l_x2;
	double l_y1, l_y2;
	double l_z1, l_z2;
}				t_rotor;


// --- CONTROLER GAIN VALS --- //
#define Kpx 40	// 50
#define Kdx 15	// 8
#define Kix 20	// 8

#define Kpy 40	// 50
#define Kdy 15	// 8
#define Kiy 20	// 8

#define Kpz 50	// 40
#define Kdz 10	// 14
#define Kiz 35	// 20

#define Kpa 45	// 80
#define Kia 5	// 10
#define Kda 5	// 8

#define Kpb 45	// 80
#define Kib 5	// 10
#define Kdb 5	// 8

#define Kpc 45	// 80
#define Kic 5	// 10 
#define Kdc 5	// 6
// ------------------------- //

// ---- Parameter ---- //
#define m 0.65
#define mp 0.3 // 0.3
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

// --- Output data file --- //
#define FD_NUM 7
#define PATH_FILE "DATA/result_line_path"
#define ERROR_FILE "DATA/result_line_error"
#define ANI_PATH_FILE "../ANIMATION_PID/DATA/path_DATA"
#define ANI_CABLE_FILE "../ANIMATION_PID/DATA/cable_DATA"
#define ANI_XROTOR_FILE "../ANIMATION_PID/DATA/xrotor_DATA"
#define ANI_YROTOR_FILE "../ANIMATION_PID/DATA/yrotor_DATA"
#define THRUST_TORQU "DATA/rotor_thrust_and_torqu"
// ------------------------ //

void	my_free(int rev_f, int size, double **ptr);
void	my_fclose(int rev_f, int size, FILE **fd);
void	set_desire_path(double t, t_desire *des);
int		func(double t, double *x, double *k, int i, t_rotor *rot, t_desire *des, t_integral *intg);

#endif
