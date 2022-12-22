/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   my_header.h                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: dmogi <mogi.t218140@gmail.com>             +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2022/12/22 18:00:15 by dmogi             #+#    #+#             */
/*   Updated: 2022/12/22 18:00:16 by dmogi            ###   ########.jp       */
/*                                                                            */
/* ************************************************************************** */

#ifndef MY_HEADER_H
# define MY_HEADER_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
#define dt 5.0e-4

// --- CONTROLER GAIN VALS --- //
#define Kpt 80
#define Kdt 60

#define Kpa 200	// 80
#define Kia 10	// 10
#define Kda 8	// 8

#define Kpb 200	// 80
#define Kib 10	// 10
#define Kdb 8	// 8

#define Kpc 80	// 80
#define Kic 10	// 10 
#define Kdc 8	// 6
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

void my_free(int size, double **ptr);
void my_fclose(int size, FILE **fd);

#endif
