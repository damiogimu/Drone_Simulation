#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_val_define.h"

#define FD_NUM 7
#define PATH_FILE "DATA/result_line_path"
#define ERROR_FILE "DATA/result_line_error"
#define ANI_PATH_FILE "../ANIMATION_PID/DATA/path_DATA"
#define ANI_CABLE_FILE "../ANIMATION_PID/DATA/cable_DATA"
#define ANI_XROTOR_FILE "../ANIMATION_PID/DATA/xrotor_DATA"
#define ANI_YROTOR_FILE "../ANIMATION_PID/DATA/yrotor_DATA"
#define THRUST_TORQU "DATA/rotor_thrust_and_torqu"

#define SCALE 1.0

static double e[4][6];
static double e_old[4][6];
static double integral[4][6];

int func(double t, double *x, double *k, int i, t_subdata *subdata)
{
	double x_dir = (cos(x[3])*sin(x[4])*cos(x[5])) + (sin(x[3])*sin(x[5]));
	double y_dir = (cos(x[3])*sin(x[4])*sin(x[5])) - (sin(x[3])*cos(x[5]));
	double z_dir = cos(x[3])*cos(x[4]);

	double xd, yd, zd, dxd, dyd, dzd, ddxd, ddyd, ddzd;
	double ad, bd, cd, dad, dbd, dcd, ddad, ddbd, ddcd;
	double ele_ad, ele_bd, dad_nume, dbd_nume, ddad_ele1, ddad_ele2, ddbd_ele1, ddbd_ele2;
	double Fx_d, Fy_d, Fz_d, U1_d, U2_d, U3_d, U4_d;

	double T1_d, T2_d, T3_d, T4_d;
	double thrust, x_mom, y_mom, z_mom;
	double coefi = (Cq/Ct)*(R/M_PI);

//---	Define Desire Path	---//

	if (t <= 2.0)
	{
		xd = 0.0;
		dxd = 0.0;
		ddxd = 0.0;

		yd = 0.0;
		dyd = 0.0;
		ddyd = 0.0;
	}
	else
	{
		t -= 2.0;
		xd = 0.05*t*t;
		dxd = 0.1*t;
		ddxd = 0.1;

		yd = 0.05*t*t;
		dyd = 0.1*t;
		ddyd = 0.1;
	}
	zd = 1.0;
	dzd = 0.0;
	ddzd = 0.0;	
	// -------- //	
	zd += 1.0e-8;  // 高度の初期位置と目標位置が一致して Fz_d=0 となるのを防ぐ	
	cd = 0.0;
	dcd = 0.0;
	ddcd = 0.0;

//----------------------------//
	subdata->des_x = xd;
	subdata->des_y = yd;
	subdata->des_z = zd;

	e[i][0] = xd - x[0];
	e[i][1] = yd - x[1];
	e[i][2] = zd - x[2];
	integral[i][0] += ((e[i][0] + e_old[i][0])/2.0) * dt;
	integral[i][1] += ((e[i][1] + e_old[i][1])/2.0) * dt;
	integral[i][2] += ((e[i][2] + e_old[i][2])/2.0) * dt;
	e_old[i][0] = e[i][0];
	e_old[i][1] = e[i][1];
	e_old[i][2] = e[i][2];

	Fx_d = Kix*integral[i][0] + Kdx*(dxd - x[8]) + Kpx*(xd - x[0]);
	Fy_d = Kiy*integral[i][1] + Kdy*(dyd - x[9]) + Kpy*(yd - x[1]);
	Fz_d = Kiz*integral[i][2] + Kdz*(dzd - x[10]) + Kpz*(zd - x[2]);

	U1_d = (Fz_d)/(cos(x[3])*cos(x[4]));
	ele_ad = (1.0/U1_d)*(Fx_d*sin(cd) - Fy_d*cos(cd));
	if (ele_ad <= -1.0)
		ele_ad = -0.999;
	else if (ele_ad >= 1.0)
		ele_ad = 0.999;
	ad = asin(ele_ad);
	ele_bd = (1.0/(U1_d*cos(ad)))*(Fx_d*cos(cd) + Fy_d*sin(cd));
	if (ele_bd <= -1.0)
		ele_bd = -0.999;
	else if (ele_bd >= 1.0)
		ele_bd = 0.999;
	bd = asin(ele_bd);

	e[i][3] = ad - x[3];
	e[i][4] = bd - x[4];
	e[i][5] = cd - x[5];
	integral[i][3] += ((e[i][3] + e_old[i][3])/2.0) * dt;
	integral[i][4] += ((e[i][4] + e_old[i][4])/2.0) * dt;
	integral[i][5] += ((e[i][5] + e_old[i][5])/2.0) * dt;
	e_old[i][3] = e[i][3];
	e_old[i][4] = e[i][4];
	e_old[i][5] = e[i][5];

	dad_nume = (Fy_d*sin(cd)*dcd + Fx_d*cos(cd)*dcd);
	dbd_nume = (Fy_d*cos(cd)*dcd - Fx_d*sin(cd)*dcd);
	ddad_ele1 = (Fy_d*sin(cd)*ddcd + Fx_d*cos(cd)*ddcd - Fx_d*sin(cd)*ddcd + Fy_d*cos(cd)*ddcd)/(U1_d*sqrt(1 - ele_ad*ele_ad));
	ddad_ele2 = ((Fx_d*sin(cd) - Fy_d*cos(cd))*(Fy_d*sin(cd)*dcd + Fx_d*cos(cd)*dcd*dcd))/(U1_d*U1_d*U1_d*pow((1-ele_ad*ele_ad), 3.0/2.0));
	ddbd_ele1 = (-Fx_d*sin(cd)*ddcd + Fy_d*cos(cd)*ddcd - Fy_d*sin(cd)*dcd*dcd - Fx_d*cos(cd)*dcd*dcd)/(U1_d*cos(cd)*sqrt(1-ele_bd*ele_bd));
	ddbd_ele2 = ((Fy_d*sin(cd) + Fx_d*cos(cd))*(Fy_d*cos(cd)*dcd-Fx_d*sin(cd)*dcd*dcd))/(U1_d*U1_d*U1_d*pow((1-ele_bd*ele_bd), 3.0/2.0)*cos(cd)*cos(cd)*cos(cd));

	dad = dad_nume/U1_d*(sqrt(1 - ele_ad*ele_ad));
	dbd = dbd_nume/(U1_d*cos(cd)*sqrt(1 - ele_bd*ele_bd));
	ddad = ddad_ele1 + ddad_ele2;
	ddbd = ddbd_ele1 + ddbd_ele2;

	U2_d = Kia*integral[i][3] + Kda*(dad - x[11]) + Kpa*(ad - x[3]);
	U3_d = Kib*integral[i][4] + Kdb*(dbd - x[12]) + Kpb*(bd - x[4]);
	U4_d = Kic*integral[i][5] + Kdc*(dcd - x[13]) + Kpc*(cd - x[5]);

	T3_d = -(U4_d/coefi + 2.0*U3_d/l - U1_d)/4.0;
	T2_d = U1_d/2.0 - U2_d/(2.0*l) - U3_d/(2.0*l) - T3_d;
	T1_d = U3_d/l + T3_d;
	T4_d = U2_d/l + T2_d;

	subdata->T1 = T1_d;
	subdata->T2 = T2_d;
	subdata->T3 = T3_d;
	subdata->T4 = T4_d;

	thrust = T1_d + T2_d + T3_d + T4_d;
	x_mom = l*(T4_d - T2_d);
	y_mom = l*(T1_d - T3_d);
	z_mom = coefi*(- T1_d + T2_d - T3_d + T4_d);

	k[0] = x[8];
	k[1] = x[9];
	k[2] = x[10];
	k[3] = x[11];
	k[4] = x[12];
	k[5] = x[13];
	k[6] = x[14];
	k[7] = x[15];
	
	k[8] = pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),3)*sin(x[6])*pow(x[14],2)-x_dir*thrust*mp*pow(cos(x[7]),2)*pow(sin(x[6]),2)+((2*g*pow(mp,2)+(2*g*m-z_dir*thrust)*mp)*pow(cos(x[7]),2)*cos(x[6])+m*mp*r*cos(x[7])*pow(x[15],2)-y_dir*thrust*mp*cos(x[7])*sin(x[7]))*sin(x[6])+x_dir*thrust*mp+x_dir*thrust*m);
	k[9] = pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),2)*sin(x[7])*pow(x[14],2)-x_dir*thrust*mp*cos(x[7])*sin(x[7])*sin(x[6])+(2*g*pow(mp,2)+(2*g*m-z_dir*thrust)*mp)*cos(x[7])*sin(x[7])*cos(x[6])+m*mp*r*sin(x[7])*pow(x[15],2)+y_dir*thrust*mp*pow(cos(x[7]),2)+y_dir*thrust*m);
	k[10] = pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),3)*cos(x[6])*pow(x[14],2)-x_dir*thrust*mp*pow(cos(x[7]),2)*cos(x[6])*sin(x[6])+(2*g*pow(mp,2)+(2*g*m-z_dir*thrust)*mp)*pow(cos(x[7]),2)*pow(cos(x[6]),2)+(m*mp*r*cos(x[7])*pow(x[15],2)-y_dir*thrust*mp*cos(x[7])*sin(x[7]))*cos(x[6])-2*g*pow(mp,2)+(z_dir*thrust-3*g*m)*mp-g*pow(m,2)+z_dir*thrust*m);

	k[11] = -pow(Ixx,-1)*pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-2)*((((Iyy-Ixx)*pow(Izz,2)+(pow(Ixx,2)-pow(Iyy,2))*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),4)+(Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*pow(cos(x[4]),2))*cos(x[3])*sin(x[3])*pow(x[13],2)+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),2)*sin(x[4])*cos(x[3])*sin(x[3])*x[11]+(((2*Iyy-Ixx)*pow(Izz,2)+(pow(Ixx,2)-2*pow(Iyy,2))*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),3)+(Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*cos(x[4]))*x[12]*pow(cos(x[3]),2)+(((Ixx-Iyy)*pow(Izz,2)+(pow(Iyy,2)-pow(Ixx,2))*Izz)*pow(cos(x[4]),3)+((pow(Ixx,2)-Ixx*Iyy)*Izz-Ixx*pow(Izz,2))*cos(x[4]))*x[12])*x[13]+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*cos(x[4])*sin(x[4])*x[12]*pow(cos(x[3]),2)+(Ixx*pow(Izz,2)+((-Ixx*Iyy)-pow(Ixx,2))*Izz)*cos(x[4])*sin(x[4])*x[12])*x[11]+((pow(Iyy,2)*Izz-Iyy*pow(Izz,2))*pow(cos(x[4]),2)*pow(x[12],2)+(Ixx*Iyy-Ixx*Izz)*y_mom*cos(x[4])*sin(x[4]))*cos(x[3])*sin(x[3])+((Ixx*Izz-Ixx*Iyy)*z_mom*sin(x[4])+(Ixx*Iyy-Ixx*Izz)*x_mom*pow(cos(x[4]),2)+(Ixx*Izz-Ixx*Iyy)*x_mom)*pow(cos(x[3]),2)-Ixx*Izz*z_mom*sin(x[4])+(Ixx-Iyy)*Izz*x_mom*pow(cos(x[4]),2)-Ixx*Izz*x_mom);
	k[12] = pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-1)*(((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*pow(cos(x[4]),2)*sin(x[4])*pow(sin(x[3]),2)+(Ixx*Izz-pow(Izz,2))*pow(cos(x[4]),2)*sin(x[4]))*pow(x[13],2)+((((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*pow(cos(x[4]),2)*pow(sin(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*pow(cos(x[4]),2))*x[11]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*x[12]*cos(x[3])*sin(x[3]))*x[13]+((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*x[12]*cos(x[3])*sin(x[3])*x[11]+(Iyy-Izz)*y_mom*cos(x[4])*pow(sin(x[3]),2)+((Izz-Iyy)*x_mom*sin(x[4])+(Izz-Iyy)*z_mom)*cos(x[3])*sin(x[3])+Izz*y_mom*cos(x[4]));
	k[13] = -pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-2)*((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*pow(cos(x[4]),2)*sin(x[4])*cos(x[3])*sin(x[3])*pow(x[13],2)+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*pow(cos(x[4]),2)*cos(x[3])*sin(x[3])*x[11]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*x[12]*pow(cos(x[3]),2)+((Ixx-Iyy)*Izz-pow(Izz,2))*cos(x[4])*sin(x[4])*x[12])*x[13]+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*x[12]*pow(cos(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*cos(x[4])*x[12])*x[11]+(Iyy-Izz)*y_mom*cos(x[4])*cos(x[3])*sin(x[3])+((Izz-Iyy)*x_mom*sin(x[4])+(Izz-Iyy)*z_mom)*pow(cos(x[3]),2)-Izz*x_mom*sin(x[4])-Izz*z_mom);

//	k[14] = 0.0;
//	k[15] = 0.0;

	k[14] = pow(m,-1)*pow(r,-1)*pow(cos(x[7]),-1)*(2*m*r*sin(x[7])*x[15]*x[14]+((-2*g*mp)-2*g*m+z_dir*thrust)*sin(x[6])-x_dir*thrust*cos(x[6]));
	k[15] = -pow(m,-1)*pow(r,-1)*(m*r*cos(x[7])*sin(x[7])*pow(x[14],2)-x_dir*thrust*sin(x[7])*sin(x[6])+(2*g*mp+2*g*m-z_dir*thrust)*sin(x[7])*cos(x[6])+y_dir*thrust*cos(x[7]));

	return (0);
}

int main(void)
{
	char files[FD_NUM][100]={PATH_FILE, ERROR_FILE, ANI_PATH_FILE, ANI_CABLE_FILE, ANI_XROTOR_FILE, ANI_YROTOR_FILE, THRUST_TORQU};
	int i, j;
	FILE *fd[FD_NUM];
	double t = 0.0;
	double x[X_NUM] = {X_POS, Y_POS, Z_POS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double c[RK4_SIZE] = {0.0, 0.5, 0.5, 1.0};
	double **state;
	double **k;
	double x_p, y_p, z_p;
	double r_x1, r_x2, r_y1, r_y2, r_z1, r_z2;
	double ddx, ddy, ddz;
	t_subdata subdata;

// --- setup --- //
	i = 0;
	while (i < FD_NUM)
	{
		if ((fd[i] = fopen(files[i], "w")) == NULL)
		{
			while (i > 0)
			{
				i--;
				fclose(fd[i]);
			}
			return (-1);
		}
		i++;
	}
	state = malloc(sizeof(double *) * RK4_SIZE);
	k = malloc(sizeof(double *) * RK4_SIZE);
	i = 0;
	while(i < RK4_SIZE)
	{
		state[i] = malloc(sizeof(double) * X_NUM);
		k[i] = malloc(sizeof(double) * X_NUM);
		i++;
	}
	i = 0;
	while(i < 4)
	{
		j = 0;
		while (j < 6)
		{
			e[i][j] = 0.0;
			e_old[i][j] = 0.0;
			integral[i][j] = 0.0;
			j++;
		}
		j = 0;
		while (j < X_NUM)
		{
			k[i][j] = 0.0;
			j++;
		}
		i++;
	}
// ------------- //
	while (t <= (TIME))
	{
		x_p = x[0] + r*cos(x[7])*sin(x[6]);
		y_p = x[1] + r*sin(x[7]);
		z_p = x[2] - r*cos(x[7])*cos(x[6]);
//-----------  output data  -----------//
		fprintf(fd[0], "%f %f %f %f %f %f %f %f %f", t, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]);
		fprintf(fd[0], " %f %f %f %f %f %f %f %f", x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
		fprintf(fd[0], " %f %f %f\n", x_p, y_p, z_p);
		if ((int)((t*1000000.0)+0.5)%50000 == 0)
		{
			fprintf(fd[3], "%f %f %f %f\n", t, x[0], x[1], x[2]);
			fprintf(fd[3], "%f %f %f %f\n", t, x_p, y_p, z_p);
			fprintf(fd[3], "\n\n");

			r_x1 = x[0]+(SCALE*l)*cos(x[5])*cos(x[5]);
			r_y1 = x[1]+(SCALE*l)*sin(x[5])*cos(x[5]);
			r_z1 = x[2] -(SCALE*l)*tan(x[4]);
			fprintf(fd[4], "%f %f %f %f\n", t, r_x1, r_y1, r_z1);
			r_x2 = x[0]-(SCALE*l)*cos(x[5])*cos(x[5]);
			r_y2 = x[1]-(SCALE*l*cos(x[5])*sin(x[5]));
			r_z2 = x[2] + (SCALE*l)*tan(x[4]);
			fprintf(fd[4], "%f %f %f %f\n", t, r_x2, r_y2, r_z2);
			fprintf(fd[4], "\n\n");

			r_x1 = x[0]+(SCALE*l*cos(x[5])*sin(x[5]));
			r_y1 = x[1]-(SCALE*l)*cos(x[5])*cos(x[5]);
			r_z1 = x[2] - (SCALE*l)*tan(x[3]);
			fprintf(fd[5], "%f %f %f %f\n", t, r_x1, r_y1, r_z1);
			r_x2 = x[0]-(SCALE*l*sin(x[5])*cos(x[5]));
			r_y2 = x[1]+(SCALE*l)*cos(x[5])*cos(x[5]);
			r_z2 = x[2] + (SCALE*l)*tan(x[3]);
			fprintf(fd[5], "%f %f %f %f\n", t, r_x2, r_y2, r_z2);
			fprintf(fd[5], "\n\n");

			fprintf(fd[2], "%f %f %f %f %f %f %f %f %f %f", t, x[0],x[1],x[2],x[3], x[4], x[5], x_p, y_p, z_p);
		}
//-------------------------------------// 
		i = 0;
		while (i < RK4_SIZE)
		{
			j = 0;
			while (j < X_NUM)
			{
				if (i == 0)
					state[i][j] = x[j];
				else
					state[i][j] = x[j] + (state[i-1][j]*dt*c[i]);
				j++;
			}
			func(t, state[i], k[i], i, &subdata);
			i++;
		}
		j = 0;
		while (j < X_NUM)
		{
			x[j] += ((k[0][j] + 2.0*k[1][j] + 2.0*k[2][j] + k[3][j])*(dt))/6.0;
			j++;
		}
		ddx = (k[0][8] + 2.0*k[1][8] + 2.0*k[2][8] + k[3][8])/6.0;
		ddy = (k[0][9] + 2.0*k[1][9] + 2.0*k[2][9] + k[3][9])/6.0;
		ddz = (k[0][10] + 2.0*k[1][10] + 2.0*k[2][10] + k[3][10])/6.0;
//-----------  output data  -----------//
		if ((int)((t*1000000.0)+0.5)%50000==0)
			fprintf(fd[2], " %f %f %f\n", subdata.des_x, subdata.des_y, subdata.des_z);
		fprintf(fd[6], "%f %f %f %f %f\n", t, subdata.T1, subdata.T1, subdata.T3, subdata.T4);
		fprintf(fd[1], "%f %f %f %f ", t, (subdata.des_x-x[0]), (subdata.des_y-x[1]), (subdata.des_z-x[2]));
		fprintf(fd[1], "%f %f %f\n", ddx, ddy, ddz);
//------------------------------------//
		t += dt;
	}
	i = 0;
	while (i < FD_NUM)
	{
		fclose(fd[i]);
		i++;
	}	
	i = 0;
	while(i < RK4_SIZE)
	{
		free(state[i]);
		free(k[i]);
		i++;
	}
	free(state);
	free(k);
	return (0);
}
