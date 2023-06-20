#include "my_header.h"

#define SCALE 1.0

int func(double t, double *x, double *k, int i, t_desire *des, t_gain *gain)
{
	double f;
	double x_dir = (cos(x[3])*sin(x[4])*cos(x[5])) + (sin(x[3])*sin(x[5]));
	double y_dir = (cos(x[3])*sin(x[4])*sin(x[5])) - (sin(x[3])*cos(x[5]));
	double z_dir = cos(x[3])*cos(x[4]);

	double ele_ad, ele_bd, dad_nume, dbd_nume, ddad_ele1, ddad_ele2, ddbd_ele1, ddbd_ele2;
	double Fx_d, Fy_d, Fz_d, U1_d, U2_d, U3_d, U4_d;

	f = sqrt(g/r)/(2.0*M_PI);
//	set_desire_path1(t, f, des);
	set_desire_path2(t, f, des);

	double O, P, Q;
	O = des->ddxd + gain->Kdt*(des->dxd-x[8]) + gain->Kpt*(des->xd-x[0]);
	P = des->ddyd + gain->Kdt*(des->dyd-x[9]) + gain->Kpt*(des->yd-x[1]);
	Q = des->ddzd + gain->Kdt*(des->dzd-x[10]) + gain->Kpt*(des->zd-x[2]);

	Fx_d = (-mp*r*pow(cos(x[7]),3)*sin(x[6])*pow(x[14],2))+mp*O*pow(cos(x[7]),2)*pow(sin(x[6]),2)+(((-mp*Q)-g*mp)*pow(cos(x[7]),2)*cos(x[6])-mp*r*cos(x[7])*pow(x[15],2)+mp*P*cos(x[7])*sin(x[7]))*sin(x[6])+m*O;
	Fy_d = (-mp*r*pow(cos(x[7]),2)*sin(x[7])*pow(x[14],2))+mp*O*cos(x[7])*sin(x[7])*sin(x[6])+((-mp*Q)-g*mp)*cos(x[7])*sin(x[7])*cos(x[6])-mp*r*sin(x[7])*pow(x[15],2)-mp*P*pow(cos(x[7]),2)+(mp+m)*P;
	Fz_d = mp*r*pow(cos(x[7]),3)*cos(x[6])*pow(x[14],2)-mp*O*pow(cos(x[7]),2)*cos(x[6])*sin(x[6])+(mp*Q+g*mp)*pow(cos(x[7]),2)*pow(cos(x[6]),2)+(mp*r*cos(x[7])*pow(x[15],2)-mp*P*cos(x[7])*sin(x[7]))*cos(x[6])+m*Q+g*m;

	U1_d = (Fz_d)/(cos(x[3])*cos(x[4]));
	ele_ad = (1.0/U1_d)*(Fx_d*sin(des->cd) - Fy_d*cos(des->cd));
	des->ad = asin(ele_ad);
	ele_bd = (1.0/(U1_d*cos(des->ad)))*(Fx_d*cos(des->cd) + Fy_d*sin(des->cd));
	des->bd = asin(ele_bd);
	des->dad = (des->ad - des->ad_old[i])/dt;
	des->dbd = (des->bd - des->bd_old[i])/dt;
	des->ddad = (des->dad - des->dad_old[i])/dt;
	des->ddbd = (des->dbd - des->dbd_old[i])/dt;
	des->ad_old[i] = des->ad;
	des->dad_old[i] = des->dad;
	des->bd_old[i] = des->bd;
	des->dbd_old[i] = des->dbd;

	U2_d = Ixx*(des->ddad+gain->Kdr*(des->dad-x[11])+gain->Kpr*(des->ad-x[3])) - Ixx*sin(x[4])*(des->ddcd+gain->Kdr*(des->dcd-x[13])+gain->Kpr*(des->cd-x[5])) + (Izz - Iyy)*cos(x[3])*sin(x[3])*cos(x[4])*cos(x[4])*x[13]*x[13] + (Iyy-Izz)*cos(x[3])*sin(x[3])*x[12]*x[12] + ((Iyy-Izz)*sin(x[3])*sin(x[3])-(Iyy-Izz)*cos(x[3])*cos(x[3])-Ixx)*cos(x[4])*x[12]*x[13];
	U3_d = (Iyy*cos(x[3])*cos(x[3]) + Izz*sin(x[3])*sin(x[3]))*(des->ddbd+gain->Kdr*(des->dbd-x[12])+gain->Kpr*(des->bd-x[4])) + (Iyy-Izz)*cos(x[3])*sin(x[3])*cos(x[4])*(des->ddcd+gain->Kdr*(des->dcd-x[13])+gain->Kpr*(des->cd-x[5])) + (Ixx - Iyy*sin(x[3])*sin(x[3]) - Izz*cos(x[3])*cos(x[3]))*sin(x[4])*cos(x[4])*x[13]*x[13] + 2*(Iyy-Izz)*sin(x[3])*cos(x[3])*x[11]*x[12] + ((Izz-Iyy)*cos(x[3])*cos(x[3]) - (Izz - Iyy)*sin(x[3])*sin(x[3]) - Ixx)*cos(x[4])*x[11]*x[13];
	U4_d = (Ixx*sin(x[4])*sin(x[4]) + Iyy*sin(x[3])*sin(x[3])*cos(x[4])*cos(x[4]) + Izz*cos(x[3])*cos(x[3])*cos(x[4])*cos(x[4]))*(des->ddcd+gain->Kdr*(des->dcd-x[13])+gain->Kpr*(des->cd-x[5])) - Ixx*sin(x[4])*(des->ddad+gain->Kdr*(des->dad-x[11])+gain->Kpr*(des->ad-x[3])) + (Iyy-Izz)*cos(x[3])*sin(x[3])*cos(x[4])*(des->ddbd+gain->Kdr*(des->dbd-x[12])+gain->Kpr*(des->bd-x[4])) + 2*(Iyy-Izz)*cos(x[3])*sin(x[3])*cos(x[4])*cos(x[4])*x[11]*x[13] + 2*(Iyy-Izz)*cos(x[4])*sin(x[4])*cos(x[3])*cos(x[3])*x[12]*x[13] + 2*(Ixx-Iyy)*cos(x[4])*sin(x[4])*x[12]*x[13] + 2*(Iyy-Izz)*cos(x[3])*cos(x[3])*cos(x[4])*x[11]*x[12] + (-Ixx-Iyy+Izz)*cos(x[4])*x[11]*x[12] + (Izz-Iyy)*sin(x[4])*cos(x[3])*sin(x[3])*x[12]*x[12];

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
	k[14] = pow(m,-1)*pow(r,-1)*pow(cos(x[7]),-1)*(2*m*r*sin(x[7])*x[15]*x[14]+((-2*g*mp)-2*g*m+z_dir*U1_d)*sin(x[6])-x_dir*U1_d*cos(x[6]));
	k[15] = -pow(m,-1)*pow(r,-1)*(m*r*cos(x[7])*sin(x[7])*pow(x[14],2)-x_dir*U1_d*sin(x[7])*sin(x[6])+(2*g*mp+2*g*m-z_dir*U1_d)*sin(x[7])*cos(x[6])+y_dir*U1_d*cos(x[7]));

	return (0);
}

int main(void)
{
	FILE *fd;
	char filename[100];
	char rm_str[103];
	int i, j;
	double t = 0.0;
	double c[RK4_SIZE] = {0.0, 0.5, 0.5, 1.0};
	double *x;
	double **state = NULL;
	double **k = NULL;
	double ddx, ddy, ddz;
	t_desire des;
	t_gain gain;

	if (setup(&state, &k, &x) == -1)
		return (0);
	init(state, k, x, &des);

	gain.Kpt = MIN_Kpt_GAIN;
	while (gain.Kpt <= MAX_Kpt_GAIN)
	{
		gain.Kdt = MIN_Kdt_GAIN;
		while (gain.Kdt <= MAX_Kdt_GAIN)
		{
			gain.Kpr = MIN_Kpr_GAIN;
			while (gain.Kpr <= MAX_Kpr_GAIN)
			{
				gain.Kdr = MIN_Kdr_GAIN;
				while (gain.Kdr <= MAX_Kdr_GAIN)
				{
					sprintf(filename, "DATA/gain%d_%d_%d_%d", gain.Kpt, gain.Kdt, gain.Kpr, gain.Kdr);
					fd = fopen(filename, "w");
					t = 0.0;
					init(state, k, x, &des);
					while (t <= TIME)
					{
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
							func(t, state[i], k[i], i, &des, &gain);
							i++;
						}
						j = 0;
						while (j < X_NUM)
						{
							x[j] += ((k[0][j] + 2.0*k[1][j] + 2.0*k[2][j] + k[3][j])*dt)/6.0;
							j++;
						}
				//-----------  output data  -----------//
						ddx = (k[0][8] + 2.0*k[1][8] + 2.0*k[2][8] + k[3][8])/6.0;
						ddy = (k[0][9] + 2.0*k[1][9] + 2.0*k[2][9] + k[3][9])/6.0;
						ddz = (k[0][10] + 2.0*k[1][10] + 2.0*k[2][10] + k[3][10])/6.0;
						fprintf(fd, "%f %f %f %f %f\n", t, ddx, ddy, des.ddxd, des.ddyd);
						if (isnan(ddx) != 0 || (4.5 <= t && 0.1 <= ddx))
						{
							sprintf(rm_str, "rm %s", filename);
							system(rm_str);
							break;
						}
				//------------------------------------//
						t += dt;
					}
					fclose(fd);
					gain.Kdr += Kdr_ITV;
				}
				gain.Kpr += Kpr_ITV;
			}
			gain.Kdt += Kdt_ITV;
		}
		gain.Kpt += Kpt_ITV;
	}

	free(x);
	my_free(0, RK4_SIZE, state);
	my_free(0, RK4_SIZE, k);
	return (0);
}
