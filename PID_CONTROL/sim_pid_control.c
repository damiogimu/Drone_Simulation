#include "my_header.h"

#define SCALE 1.0

int func(double t, double *x, double *k, int i, t_rotor *rot, t_desire *des, t_integral *intg)
{
	double f;
	double x_dir = (cos(x[3])*sin(x[4])*cos(x[5])) + (sin(x[3])*sin(x[5]));
	double y_dir = (cos(x[3])*sin(x[4])*sin(x[5])) - (sin(x[3])*cos(x[5]));
	double z_dir = cos(x[3])*cos(x[4]);

	double ele_ad, ele_bd;
	double Fx_d, Fy_d, Fz_d, U1_d, U2_d, U3_d, U4_d;

	f = sqrt(g/r)/(2.0*M_PI);
//	any_traj(t, des);
//	noncontrol_traj1(t, des);
	controled_traj1(t, f, des);
//	noncontrol_traj2(t, des);
//	controled_traj2(t, f, des);

	intg->e[i][0] = des->xd - x[0];
	intg->e[i][1] = des->yd - x[1];
	intg->e[i][2] = des->zd - x[2];
	intg->area[i][0] += ((intg->e[i][0] + intg->e_old[i][0])/2.0)*dt;
	intg->area[i][1] += ((intg->e[i][1] + intg->e_old[i][1])/2.0)*dt;
	intg->area[i][2] += ((intg->e[i][2] + intg->e_old[i][2])/2.0)*dt;
	intg->e_old[i][0] = intg->e[i][0];
	intg->e_old[i][1] = intg->e[i][1];
	intg->e_old[i][2] = intg->e[i][2];

	Fx_d = Kit*intg->area[i][0] + Kdt*(des->dxd - x[8]) + Kpt*(des->xd - x[0]);
	Fy_d = Kit*intg->area[i][1] + Kdt*(des->dyd - x[9]) + Kpt*(des->yd - x[1]);
	Fz_d = Kiz*intg->area[i][2] + Kdz*(des->dzd - x[10]) + Kpz*(des->zd - x[2]);

	U1_d = (Fz_d)/(cos(x[3])*cos(x[4]));
	ele_ad = (1.0/U1_d)*(Fx_d*sin(des->cd) - Fy_d*cos(des->cd));
	if (ele_ad <= -1.0)
		ele_ad = -0.999;
	else if (ele_ad >= 1.0)
		ele_ad = 0.999;
	des->ad = asin(ele_ad);
	ele_bd = (1.0/(U1_d*cos(des->ad)))*(Fx_d*cos(des->cd) + Fy_d*sin(des->cd));
	if (ele_bd <= -1.0)
		ele_bd = -0.999;
	else if (ele_bd >= 1.0)
		ele_bd = 0.999;
	des->bd = asin(ele_bd);
	des->dad = (des->ad - des->ad_old[i])/dt;
	des->dbd = (des->bd - des->bd_old[i])/dt;
	des->ddad = (des->dad - des->dad_old[i])/dt;
	des->ddbd = (des->dbd - des->dbd_old[i])/dt;
	des->ad_old[i] = des->ad;
	des->bd_old[i] = des->bd;
	des->dad_old[i] = des->dad;
	des->dbd_old[i] = des->dbd;

	intg->e[i][3] = des->ad - x[3];
	intg->e[i][4] = des->bd - x[4];
	intg->e[i][5] = des->cd - x[5];
	intg->area[i][3] += ((intg->e[i][3] + intg->e_old[i][3])/2.0) * dt;
	intg->area[i][4] += ((intg->e[i][4] + intg->e_old[i][4])/2.0) * dt;
	intg->area[i][5] += ((intg->e[i][5] + intg->e_old[i][5])/2.0) * dt;
	intg->e_old[i][3] = intg->e[i][3];
	intg->e_old[i][4] = intg->e[i][4];
	intg->e_old[i][5] = intg->e[i][5];
	
	U2_d = Kir*intg->area[i][3] + Kdr*(des->dad - x[11]) + Kpr*(des->ad - x[3]);
	U3_d = Kir*intg->area[i][4] + Kdr*(des->dbd - x[12]) + Kpr*(des->bd - x[4]);
	U4_d = Kir*intg->area[i][5] + Kdr*(des->dcd - x[13]) + Kpr*(des->cd - x[5]);

	rot->T3 = -(U4_d/rot->coefi + 2.0*U3_d/l - U1_d)/4.0;
	rot->T2 = U1_d/2.0 - U2_d/(2.0*l) - U3_d/(2.0*l) - rot->T3;
	rot->T1 = U3_d/l + rot->T3;
	rot->T4 = U2_d/l + rot->T2;

	k[0] = x[8];
	k[1] = x[9];
	k[2] = x[10];
	k[3] = x[11];
	k[4] = x[12];
	k[5] = x[13];
	k[6] = x[14];
	k[7] = x[15];

	k[8] = pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),3)*sin(x[6])*pow(x[14],2)-x_dir*U1_d*mp*pow(cos(x[7]),2)*pow(sin(x[6]),2)+((2*g*pow(mp,2)+(2*g*m-z_dir*U1_d)*mp)*pow(cos(x[7]),2)*cos(x[6])+m*mp*r*cos(x[7])*pow(x[15],2)-y_dir*U1_d*mp*cos(x[7])*sin(x[7]))*sin(x[6])+x_dir*U1_d*mp+x_dir*U1_d*m);
	k[9] = pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),2)*sin(x[7])*pow(x[14],2)-x_dir*U1_d*mp*cos(x[7])*sin(x[7])*sin(x[6])+(2*g*pow(mp,2)+(2*g*m-z_dir*U1_d)*mp)*cos(x[7])*sin(x[7])*cos(x[6])+m*mp*r*sin(x[7])*pow(x[15],2)+y_dir*U1_d*mp*pow(cos(x[7]),2)+y_dir*U1_d*m);
	k[10] = pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),3)*cos(x[6])*pow(x[14],2)-x_dir*U1_d*mp*pow(cos(x[7]),2)*cos(x[6])*sin(x[6])+(2*g*pow(mp,2)+(2*g*m-z_dir*U1_d)*mp)*pow(cos(x[7]),2)*pow(cos(x[6]),2)+(m*mp*r*cos(x[7])*pow(x[15],2)-y_dir*U1_d*mp*cos(x[7])*sin(x[7]))*cos(x[6])-2*g*pow(mp,2)+(z_dir*U1_d-3*g*m)*mp-g*pow(m,2)+z_dir*U1_d*m);
	k[10] = pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),3)*cos(x[6])*pow(x[14],2)-x_dir*U1_d*mp*pow(cos(x[7]),2)*cos(x[6])*sin(x[6])+(2*g*pow(mp,2)+(2*g*m-z_dir*U1_d)*mp)*pow(cos(x[7]),2)*pow(cos(x[6]),2)+(m*mp*r*cos(x[7])*pow(x[15],2)-y_dir*U1_d*mp*cos(x[7])*sin(x[7]))*cos(x[6])-2*g*pow(mp,2)+(z_dir*U1_d-3*g*m)*mp-g*pow(m,2)+z_dir*U1_d*m);

	k[11] = -pow(Ixx,-1)*pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-2)*((((Iyy-Ixx)*pow(Izz,2)+(pow(Ixx,2)-pow(Iyy,2))*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),4)+(Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*pow(cos(x[4]),2))*cos(x[3])*sin(x[3])*pow(x[13],2)+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),2)*sin(x[4])*cos(x[3])*sin(x[3])*x[11]+(((2*Iyy-Ixx)*pow(Izz,2)+(pow(Ixx,2)-2*pow(Iyy,2))*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),3)+(Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*cos(x[4]))*x[12]*pow(cos(x[3]),2)+(((Ixx-Iyy)*pow(Izz,2)+(pow(Iyy,2)-pow(Ixx,2))*Izz)*pow(cos(x[4]),3)+((pow(Ixx,2)-Ixx*Iyy)*Izz-Ixx*pow(Izz,2))*cos(x[4]))*x[12])*x[13]+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*cos(x[4])*sin(x[4])*x[12]*pow(cos(x[3]),2)+(Ixx*pow(Izz,2)+((-Ixx*Iyy)-pow(Ixx,2))*Izz)*cos(x[4])*sin(x[4])*x[12])*x[11]+((pow(Iyy,2)*Izz-Iyy*pow(Izz,2))*pow(cos(x[4]),2)*pow(x[12],2)+(Ixx*Iyy-Ixx*Izz)*U3_d*cos(x[4])*sin(x[4]))*cos(x[3])*sin(x[3])+((Ixx*Izz-Ixx*Iyy)*U4_d*sin(x[4])+(Ixx*Iyy-Ixx*Izz)*U2_d*pow(cos(x[4]),2)+(Ixx*Izz-Ixx*Iyy)*U2_d)*pow(cos(x[3]),2)-Ixx*Izz*U4_d*sin(x[4])+(Ixx-Iyy)*Izz*U2_d*pow(cos(x[4]),2)-Ixx*Izz*U2_d);
	k[12] = pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-1)*(((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*pow(cos(x[4]),2)*sin(x[4])*pow(sin(x[3]),2)+(Ixx*Izz-pow(Izz,2))*pow(cos(x[4]),2)*sin(x[4]))*pow(x[13],2)+((((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*pow(cos(x[4]),2)*pow(sin(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*pow(cos(x[4]),2))*x[11]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*x[12]*cos(x[3])*sin(x[3]))*x[13]+((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*x[12]*cos(x[3])*sin(x[3])*x[11]+(Iyy-Izz)*U3_d*cos(x[4])*pow(sin(x[3]),2)+((Izz-Iyy)*U2_d*sin(x[4])+(Izz-Iyy)*U4_d)*cos(x[3])*sin(x[3])+Izz*U3_d*cos(x[4]));
	k[13] = -pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-2)*((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*pow(cos(x[4]),2)*sin(x[4])*cos(x[3])*sin(x[3])*pow(x[13],2)+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*pow(cos(x[4]),2)*cos(x[3])*sin(x[3])*x[11]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*x[12]*pow(cos(x[3]),2)+((Ixx-Iyy)*Izz-pow(Izz,2))*cos(x[4])*sin(x[4])*x[12])*x[13]+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*x[12]*pow(cos(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*cos(x[4])*x[12])*x[11]+(Iyy-Izz)*U3_d*cos(x[4])*cos(x[3])*sin(x[3])+((Izz-Iyy)*U2_d*sin(x[4])+(Izz-Iyy)*U4_d)*pow(cos(x[3]),2)-Izz*U2_d*sin(x[4])-Izz*U4_d);

	k[14] = pow(r,-1)*pow(cos(x[7]),-1)*(2*r*sin(x[7])*x[15]*x[14]+((-k[10])-g)*sin(x[6])-k[8]*cos(x[6]));
	k[15] = -pow(r,-1)*(r*cos(x[7])*sin(x[7])*pow(x[14],2)-k[8]*sin(x[7])*sin(x[6])+(k[10]+g)*sin(x[7])*cos(x[6])+k[9]*cos(x[7]));

	return (0);
}

int main(void)
{
	FILE *fd[FD_NUM];
	int i, j;
	double t;
	double x[X_NUM] = {INIT_X, INIT_Y, INIT_Z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double c[RK4_SIZE] = {0.0, 0.5, 0.5, 1.0};
	double **state = NULL;
	double **k = NULL;
	double x_p, y_p, z_p;
	double ddx, ddy, ddz;
	t_desire des;
	t_rotor rot;
	t_integral intg;

	if (setup(&state, &k, fd) == -1)
		return (0);
	init(k, &intg, &des, &rot);

	t = 0.0;
	while (t <= TIME)
	{
//----------  output data  -----------//
		x_p = x[0] + r*cos(x[7])*sin(x[6]);
		y_p = x[1] + r*sin(x[7]);
		z_p = x[2] - r*cos(x[7])*cos(x[6]);
		fprintf(fd[0], "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", \
						t,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12], \
						x[13],x[14],x[15],x_p,y_p,z_p);
		if ((int)((t*1000000.0)+0.5)%50000==0)
		{
			fprintf(fd[3], "%f %f %f %f\n", t, x[0], x[1], x[2]);
			fprintf(fd[3], "%f %f %f %f\n", t, x_p, y_p, z_p);
			fprintf(fd[3], "\n\n");
			rot.l_x1 = x[0]+(SCALE*l)*cos(x[5])*cos(x[5]);
			rot.l_y1 = x[1]+(SCALE*l)*sin(x[5])*cos(x[5]);
			rot.l_z1 = x[2] -(SCALE*l)*tan(x[4]);
			fprintf(fd[4], "%f %f %f %f\n", t, rot.l_x1, rot.l_y1, rot.l_z1);
			rot.l_x2 = x[0]-(SCALE*l)*cos(x[5])*cos(x[5]);
			rot.l_y2 = x[1]-(SCALE*l*cos(x[5])*sin(x[5]));
			rot.l_z2 = x[2] + (SCALE*l)*tan(x[4]);
			fprintf(fd[4], "%f %f %f %f\n", t, rot.l_x2, rot.l_y2, rot.l_z2);
			fprintf(fd[4], "\n\n");
			rot.l_x1 = x[0]+(SCALE*l*cos(x[5])*sin(x[5]));
			rot.l_y1 = x[1]-(SCALE*l)*cos(x[5])*cos(x[5]);
			rot.l_z1 = x[2] - (SCALE*l)*tan(x[3]);
			fprintf(fd[5], "%f %f %f %f\n", t, rot.l_x1, rot.l_y1, rot.l_z1);
			rot.l_x2 = x[0]-(SCALE*l*sin(x[5])*cos(x[5]));
			rot.l_y2 = x[1]+(SCALE*l)*cos(x[5])*cos(x[5]);
			rot.l_z2 = x[2] + (SCALE*l)*tan(x[3]);
			fprintf(fd[5], "%f %f %f %f\n", t, rot.l_x2, rot.l_y2, rot.l_z2);
			fprintf(fd[5], "\n\n");
			fprintf(fd[2], "%f %f %f %f %f %f %f %f %f %f", t,x[0],x[1],x[2],x[3],x[4],x[5],x_p,y_p,z_p);
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
			func(t, state[i], k[i], i, &rot, &des, &intg);
			i++;
		}
		j = 0;
		while (j < X_NUM)
		{
			x[j] += ((k[0][j] + 2.0*k[1][j] + 2.0*k[2][j] + k[3][j])*(dt))/6.0;
			j++;
		}
//-----------  output data  -----------//
		ddx = (k[0][8] + 2.0*k[1][8] + 2.0*k[2][8] + k[3][8])/6.0;
		ddy = (k[0][9] + 2.0*k[1][9] + 2.0*k[2][9] + k[3][9])/6.0;
		ddz = (k[0][10] + 2.0*k[1][10] + 2.0*k[2][10] + k[3][10])/6.0;
		fprintf(fd[7], "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, des.ddxd, des.ddyd, des.ddzd, ddx, ddy, ddz, des.dxd, des.dyd, des.dzd, des.xd, des.yd, des.zd, des.ad, des.bd, des.cd);
		if ((int)((t*1000000.0)+0.5)%50000==0)
			fprintf(fd[2], " %f %f %f\n", des.xd, des.yd, des.zd);
		fprintf(fd[6], "%f %f %f %f %f\n", t,rot.T1, rot.T2, rot.T3, rot.T4);
		fprintf(fd[1], "%f %f %f %f ", t, (des.xd-x[0]), (des.yd-x[1]), (des.zd-x[2]));
//------------------------------------//
		t += dt;
	}
	my_fclose(0, FD_NUM, fd);
	my_free(0, RK4_SIZE, state);
	my_free(0, RK4_SIZE, k);

	return (0);
}
