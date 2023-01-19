#include "my_header.h"

#define SCALE 1.0

int func(double t, double *x, double *k, int i, t_desire *des, t_integral *intg, t_gain *gain)
{
	double x_dir = (cos(x[3])*sin(x[4])*cos(x[5])) + (sin(x[3])*sin(x[5]));
	double y_dir = (cos(x[3])*sin(x[4])*sin(x[5])) - (sin(x[3])*cos(x[5]));
	double z_dir = cos(x[3])*cos(x[4]);

	double ele_ad, ele_bd, dad_nume, dbd_nume, ddad_ele1, ddad_ele2, ddbd_ele1, ddbd_ele2;
	double Fx_d, Fy_d, Fz_d, U1_d, U2_d, U3_d, U4_d;

	vertical_traj(t, des);

	intg->e[i][0] = des->xd - x[0];
	intg->e[i][1] = des->yd - x[1];
	intg->e[i][2] = des->zd - x[2];
	intg->area[i][0] += ((intg->e[i][0] + intg->e_old[i][0])/2.0) * dt;
	intg->area[i][1] += ((intg->e[i][1] + intg->e_old[i][1])/2.0) * dt;
	intg->area[i][2] += ((intg->e[i][2] + intg->e_old[i][2])/2.0) * dt;
	intg->e_old[i][0] = intg->e[i][0];
	intg->e_old[i][1] = intg->e[i][1];
	intg->e_old[i][2] = intg->e[i][2];
	
	Fx_d = 0.0;
	Fy_d = 0.0;
	Fz_d = gain->It*intg->area[i][2] + gain->Dt*(des->dzd - x[10]) + gain->Pt*(des->zd - x[2]);

	U1_d = (Fz_d)/(cos(x[3])*cos(x[4]));
	ele_ad = (1.0/U1_d)*(Fx_d*sin(des->cd) - Fy_d*cos(des->cd));
	des->ad = asin(ele_ad);
	ele_bd = (1.0/(U1_d*cos(des->ad)))*(Fx_d*cos(des->cd) + Fy_d*sin(des->cd));
	des->bd = asin(ele_bd);

	intg->e[i][3] = des->ad - x[3];
	intg->e[i][4] = des->bd - x[4];
	intg->e[i][5] = des->cd - x[5];
	intg->area[i][3] += ((intg->e[i][3] + intg->e_old[i][3])/2.0) * dt;
	intg->area[i][4] += ((intg->e[i][4] + intg->e_old[i][4])/2.0) * dt;
	intg->area[i][5] += ((intg->e[i][5] + intg->e_old[i][5])/2.0) * dt;
	intg->e_old[i][3] = intg->e[i][3];
	intg->e_old[i][4] = intg->e[i][4];
	intg->e_old[i][5] = intg->e[i][5];

	dad_nume = (Fy_d*sin(des->cd)*des->dcd + Fx_d*cos(des->cd)*des->dcd);
	dbd_nume = (Fy_d*cos(des->cd)*des->dcd - Fx_d*sin(des->cd)*des->dcd);
	ddad_ele1 = (Fy_d*sin(des->cd)*des->ddcd + Fx_d*cos(des->cd)*des->ddcd - Fx_d*sin(des->cd)*des->ddcd + Fy_d*cos(des->cd)*des->ddcd)/(U1_d*sqrt(1 - ele_ad*ele_ad));
	ddad_ele2 = ((Fx_d*sin(des->cd) - Fy_d*cos(des->cd))*(Fy_d*sin(des->cd)*des->dcd + Fx_d*cos(des->cd)*des->dcd*des->dcd))/(U1_d*U1_d*U1_d*pow((1-ele_ad*ele_ad), 3.0/2.0));
	ddbd_ele1 = (-Fx_d*sin(des->cd)*des->ddcd + Fy_d*cos(des->cd)*des->ddcd - Fy_d*sin(des->cd)*des->dcd*des->dcd - Fx_d*cos(des->cd)*des->dcd*des->dcd)/(U1_d*cos(des->cd)*sqrt(1-ele_bd*ele_bd));
	ddbd_ele2 = ((Fy_d*sin(des->cd) + Fx_d*cos(des->cd))*(Fy_d*cos(des->cd)*des->dcd-Fx_d*sin(des->cd)*des->dcd*des->dcd))/(U1_d*U1_d*U1_d*pow((1-ele_bd*ele_bd), 3.0/2.0)*cos(des->cd)*cos(des->cd)*cos(des->cd));

	des->dad = dad_nume/U1_d*(sqrt(1 - ele_ad*ele_ad));
	des->dbd = dbd_nume/(U1_d*cos(des->cd)*sqrt(1 - ele_bd*ele_bd));
	des->ddad = ddad_ele1 + ddad_ele2;
	des->ddbd = ddbd_ele1 + ddbd_ele2;

	U2_d = 0.0;
	U3_d = 0.0;
	U4_d = 0.0;

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
	t_integral intg;
	t_gain gain;

	if (setup(&state, &k, &x) == -1)
		return (0);
	init(state, k, x, &intg, &des);

	gain.Pt = MIN_Pz_GAIN;
	while (gain.Pt <= MAX_Pz_GAIN)
	{
		gain.Dt = MIN_Dz_GAIN;
		while (gain.Dt <= MAX_Dz_GAIN)
		{
			gain.It = MIN_Iz_GAIN;
			while (gain.It <= MAX_Iz_GAIN)
			{
				sprintf(filename, "Z_DATA/gain_%d_%d_%d", gain.Pt, gain.Dt, gain.It);
				fd = fopen(filename, "w");
				t = 0.0;
				init(state, k, x, &intg, &des);
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
						func(t, state[i], k[i], i, &des, &intg, &gain);
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
					fprintf(fd, "%f %f %f\n", t, x[2], des.zd);
					if (isnan(ddx) != 0 || (3.0 <= t && 0.1 <= ddx))
					{
						sprintf(rm_str, "rm %s", filename);
						system(rm_str);
						break;
					}
			//------------------------------------//
					t += dt;
				}
				fclose(fd);
				gain.It += Iz_ITV;
			}
			gain.Dt += Dz_ITV;
		}
		gain.Pt += Pz_ITV;
	}

	free(x);
	my_free(0, RK4_SIZE, state);
	my_free(0, RK4_SIZE, k);
	return (0);
}
