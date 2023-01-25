#include "my_header.h"

int func(double t, double *x, double *k, int i, t_rotor *rot, t_desire *des)
{
	double f;
	double x_dir = (cos(x[3])*sin(x[4])*cos(x[5])) + (sin(x[3])*sin(x[5]));
	double y_dir = (cos(x[3])*sin(x[4])*sin(x[5])) - (sin(x[3])*cos(x[5]));
	double z_dir = cos(x[3])*cos(x[4]);
	double ele_ad, ele_bd;
	double Fx_d, Fy_d, Fz_d, U1_d, U2_d, U3_d, U4_d;
	double ctm_x, ctm_y, ctm_z;

	f = sqrt(g/r)/(2.0*M_PI);
//	any_traj(t, des);
//	noncontrol_traj1(t, des);
//	controled_traj1(t, f, des);
//	noncontrol_traj2(t, des);
	controled_traj2(t, f, des);

	ctm_x = des->ddxd + Kdt*(des->dxd-x[8]) + Kpt*(des->xd-x[0]);
	ctm_y = des->ddyd + Kdt*(des->dyd-x[9]) + Kpt*(des->yd-x[1]);
	ctm_z = des->ddzd + Kdt*(des->dzd-x[10]) + Kpt*(des->zd-x[2]);

	Fx_d = (-mp*r*pow(cos(x[7]),3)*sin(x[6])*pow(x[14],2))+mp*ctm_x*pow(cos(x[7]),2)*pow(sin(x[6]),2)+(((-mp*ctm_z)-g*mp)*pow(cos(x[7]),2)*cos(x[6])-mp*r*cos(x[7])*pow(x[15],2)+mp*ctm_y*cos(x[7])*sin(x[7]))*sin(x[6])+m*ctm_x;
	Fy_d = (-mp*r*pow(cos(x[7]),2)*sin(x[7])*pow(x[14],2))+mp*ctm_x*cos(x[7])*sin(x[7])*sin(x[6])+((-mp*ctm_z)-g*mp)*cos(x[7])*sin(x[7])*cos(x[6])-mp*r*sin(x[7])*pow(x[15],2)-mp*ctm_y*pow(cos(x[7]),2)+(mp+m)*ctm_y;
	Fz_d = mp*r*pow(cos(x[7]),3)*cos(x[6])*pow(x[14],2)-mp*ctm_x*pow(cos(x[7]),2)*cos(x[6])*sin(x[6])+(mp*ctm_z+g*mp)*pow(cos(x[7]),2)*pow(cos(x[6]),2)+(mp*r*cos(x[7])*pow(x[15],2)-mp*ctm_y*cos(x[7])*sin(x[7]))*cos(x[6])+m*ctm_z+g*m;

	U1_d = (Fz_d)/(cos(x[3])*cos(x[4]));
	ele_ad = (1.0/U1_d)*(Fx_d*sin(des->cd) - Fy_d*cos(des->cd));
	des->ad = asin(ele_ad);
	ele_bd = (1.0/(U1_d*cos(des->ad)))*(Fx_d*cos(des->cd) + Fy_d*sin(des->cd));
	des->bd = asin(ele_bd);

	des->dad = (des->ad - des->ad_old[i])/dt;
	des->ddad = (des->dad - des->dad_old[i])/dt;
	des->dbd = (des->bd - des->bd_old[i])/dt;
	des->ddbd = (des->dbd - des->dbd_old[i])/dt;
	des->ad_old[i] = des->ad;
	des->dad_old[i] = des->dad;
	des->bd_old[i] = des->bd;
	des->dbd_old[i] = des->dbd;

	U2_d = Ixx*(des->ddad+Kdr*(des->dad-x[11])+Kpr*(des->ad-x[3])) - Ixx*sin(x[4])*(des->ddcd+Kdr*(des->dcd-x[13])+Kpr*(des->cd-x[5])) + (Izz - Iyy)*cos(x[3])*sin(x[3])*cos(x[4])*cos(x[4])*x[13]*x[13] + (Iyy-Izz)*cos(x[3])*sin(x[3])*x[12]*x[12] + ((Iyy-Izz)*sin(x[3])*sin(x[3])-(Iyy-Izz)*cos(x[3])*cos(x[3])-Ixx)*cos(x[4])*x[12]*x[13];
	U3_d = (Iyy*cos(x[3])*cos(x[3]) + Izz*sin(x[3])*sin(x[3]))*(des->ddbd+Kdr*(des->dbd-x[12])+Kpr*(des->bd-x[4])) + (Iyy-Izz)*cos(x[3])*sin(x[3])*cos(x[4])*(des->ddcd+Kdr*(des->dcd-x[13])+Kpr*(des->cd-x[5])) + (Ixx - Iyy*sin(x[3])*sin(x[3]) - Izz*cos(x[3])*cos(x[3]))*sin(x[4])*cos(x[4])*x[13]*x[13] + 2*(Iyy-Izz)*sin(x[3])*cos(x[3])*x[11]*x[12] + ((Izz-Iyy)*cos(x[3])*cos(x[3]) - (Izz - Iyy)*sin(x[3])*sin(x[3]) - Ixx)*cos(x[4])*x[11]*x[13];
	U4_d = (Ixx*sin(x[4])*sin(x[4]) + Iyy*sin(x[3])*sin(x[3])*cos(x[4])*cos(x[4]) + Izz*cos(x[3])*cos(x[3])*cos(x[4])*cos(x[4]))*(des->ddcd+Kdr*(des->dcd-x[13])+Kpr*(des->cd-x[5])) - Ixx*sin(x[4])*(des->ddad+Kdr*(des->dad-x[11])+Kpr*(des->ad-x[3])) + (Iyy-Izz)*cos(x[3])*sin(x[3])*cos(x[4])*(des->ddbd+Kdr*(des->dbd-x[12])+Kpr*(des->bd-x[4])) + 2*(Iyy-Izz)*cos(x[3])*sin(x[3])*cos(x[4])*cos(x[4])*x[11]*x[13] + 2*(Iyy-Izz)*cos(x[4])*sin(x[4])*cos(x[3])*cos(x[3])*x[12]*x[13] + 2*(Ixx-Iyy)*cos(x[4])*sin(x[4])*x[12]*x[13] + 2*(Iyy-Izz)*cos(x[3])*cos(x[3])*cos(x[4])*x[11]*x[12] + (-Ixx-Iyy+Izz)*cos(x[4])*x[11]*x[12] + (Izz-Iyy)*sin(x[4])*cos(x[3])*sin(x[3])*x[12]*x[12];

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
	t_desire des;
	t_rotor rot;

	if (setup(&state, &k, fd) == -1)
		return (0);
	init(k, &rot, &des);

	t = 0.0;
	while (t <= TIME)
	{
		output_data1(t, x, fd, &rot);
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
			func(t, state[i], k[i], i, &rot, &des);
			i++;
		}
		j = 0;
		while (j < X_NUM)
		{
			x[j] += ((k[0][j] + 2.0*k[1][j] + 2.0*k[2][j] + k[3][j])*(dt))/6.0;
			j++;
		}
		output_data2(t, x, k, fd, &rot, &des);
		t += dt;
	}
	my_fclose(0, FD_NUM, fd);
	my_free(0, RK4_SIZE, state);
	my_free(0, RK4_SIZE, k);

	return (0);
}
