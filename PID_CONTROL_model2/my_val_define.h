#ifndef MY_VAL_DEFINE_H
# define MY_VAL_DEFINE_H

typedef struct	s_subdata 
{
	double des_x;
	double des_y;
	double des_z;
	double T1;
	double T2;
	double T3;
	double T4;
	double thrust;
}				t_subdata;

#ifndef X_POS
# define X_POS 0.0
#endif
#ifndef Y_POS
# define Y_POS 0.0
#endif
#ifndef Z_POS
# define Z_POS 0.0
#endif
#ifndef TIME
# define TIME 10.0
#endif

#define RK4_SIZE 4
#define X_NUM 16
#define dt 1.0e-3

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

// ---- FALGS ---- //
#ifndef TEST
# define TEST 0
#endif
#ifndef ERROR
# define ERROR 0
#endif
#ifndef ANIME_F
# define ANIME_F 0
#endif
#ifndef CABLE
# define CABLE 0
#endif
#ifndef X_ROTOR
# define X_ROTOR 0
#endif
#ifndef Y_ROTOR
# define Y_ROTOR 0
#endif
// -------------- //

#endif
