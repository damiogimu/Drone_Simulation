#include "my_header.h"

void func(double t, double *x, double *k)
{
	double U1, U2, U3, U4;
	double x_dir = (cos(x[3])*sin(x[4])*cos(x[5])) + (sin(x[3])*sin(x[5]));
	double y_dir = (cos(x[3])*sin(x[4])*sin(x[5])) - (sin(x[3])*cos(x[5]));
	double z_dir = cos(x[3])*cos(x[4]);
/*
	double om1=239.0619045;
	double om2=239.0619045;
	double om3=239.0619045;
	double om4=239.0619045;
	U1 = (Ct*4*p*R*R*R*R)*(om1*om1 + om2*om2 + om3*om3 + om4*om4)/(M_PI*M_PI);
	U2 = l*(Ct*4*p*R*R*R*R)*(om4*om4 - om2*om2)/(M_PI*M_PI);
	U3 = l*(Ct*4*p*R*R*R*R)*(om3*om3 - om1*om1)/(M_PI*M_PI);
	U4 = (Cq*4*p*R*R*R*R*R)*(-(om1*om1) + (om2*om2) - (om3*om3) + (om4*om4))/(M_PI*M_PI*M_PI);
*/
	U1 = 0.2*(m+mp) + (m+mp)*g;
	U2 = 0.0;
	U3 = 0.0;
	U4 = 0.0;
	
	k[0] = x[8];
	k[1] = x[9];
	k[2] = x[10];
	k[3] = x[11];
	k[4] = x[12];
	k[5] = x[13];
	k[6] = x[14];
	k[7] = x[15];

	k[8] = pow(mp+m,-1)*((((3*mp*r*cos(x[7])*x[15]+m*x[9])*sin(x[4])*sin(x[3])+(3*mp*r*cos(x[7])*sin(x[6])*x[14]+3*mp*r*sin(x[7])*x[15]*cos(x[6])+m*x[10])*sin(x[4])*cos(x[3])+(3*mp*r*cos(x[7])*cos(x[6])*x[14]-3*mp*r*sin(x[7])*x[15]*sin(x[6])+m*x[8])*cos(x[4]))*sin(x[5])+(((-3*mp*r*cos(x[7])*sin(x[6])*x[14])-3*mp*r*sin(x[7])*x[15]*cos(x[6])-m*x[10])*sin(x[3])+(3*mp*r*cos(x[7])*x[15]+m*x[9])*cos(x[3]))*cos(x[5]))*x[13]+((((-3*mp*r*cos(x[7])*x[15])-m*x[9])*sin(x[3])+((-3*mp*r*cos(x[7])*sin(x[6])*x[14])-3*mp*r*sin(x[7])*x[15]*cos(x[6])-m*x[10])*cos(x[3]))*x[11]+((-mp*r*cos(x[7])*sin(x[6])*k[14])-mp*r*cos(x[7])*cos(x[6])*pow(x[14],2)+2*mp*r*sin(x[7])*x[15]*sin(x[6])*x[14]+((-mp*r*sin(x[7])*k[15])-mp*r*cos(x[7])*pow(x[15],2))*cos(x[6]))*sin(x[3])+(mp*r*cos(x[7])*k[15]-mp*r*sin(x[7])*pow(x[15],2))*cos(x[3]))*sin(x[5])+(((3*mp*r*cos(x[7])*sin(x[6])*x[14]+3*mp*r*sin(x[7])*x[15]*cos(x[6])+m*x[10])*sin(x[4])*sin(x[3])+((-3*mp*r*cos(x[7])*x[15])-m*x[9])*sin(x[4])*cos(x[3]))*x[11]+(((-3*mp*r*cos(x[7])*x[15])-m*x[9])*cos(x[4])*x[12]+(mp*r*sin(x[7])*pow(x[15],2)-mp*r*cos(x[7])*k[15])*sin(x[4]))*sin(x[3])+(((-3*mp*r*cos(x[7])*sin(x[6])*x[14])-3*mp*r*sin(x[7])*x[15]*cos(x[6])-m*x[10])*cos(x[4])*x[12]+((-mp*r*cos(x[7])*sin(x[6])*k[14])-mp*r*cos(x[7])*cos(x[6])*pow(x[14],2)+2*mp*r*sin(x[7])*x[15]*sin(x[6])*x[14]+((-mp*r*sin(x[7])*k[15])-mp*r*cos(x[7])*pow(x[15],2))*cos(x[6]))*sin(x[4]))*cos(x[3])+(3*mp*r*cos(x[7])*cos(x[6])*x[14]-3*mp*r*sin(x[7])*x[15]*sin(x[6])+m*x[8])*sin(x[4])*x[12]+((-mp*r*cos(x[7])*cos(x[6])*k[14])+mp*r*cos(x[7])*sin(x[6])*pow(x[14],2)+2*mp*r*sin(x[7])*x[15]*cos(x[6])*x[14]+(mp*r*sin(x[7])*k[15]+mp*r*cos(x[7])*pow(x[15],2))*sin(x[6]))*cos(x[4]))*cos(x[5])+x_dir*U1);
	k[9] = -pow(mp+m,-1)*((((3*mp*r*cos(x[7])*sin(x[6])*x[14]+3*mp*r*sin(x[7])*x[15]*cos(x[6])+m*x[10])*sin(x[3])+((-3*mp*r*cos(x[7])*x[15])-m*x[9])*cos(x[3]))*sin(x[5])+((3*mp*r*cos(x[7])*x[15]+m*x[9])*sin(x[4])*sin(x[3])+(3*mp*r*cos(x[7])*sin(x[6])*x[14]+3*mp*r*sin(x[7])*x[15]*cos(x[6])+m*x[10])*sin(x[4])*cos(x[3])+(3*mp*r*cos(x[7])*cos(x[6])*x[14]-3*mp*r*sin(x[7])*x[15]*sin(x[6])+m*x[8])*cos(x[4]))*cos(x[5]))*x[13]+((((-3*mp*r*cos(x[7])*sin(x[6])*x[14])-3*mp*r*sin(x[7])*x[15]*cos(x[6])-m*x[10])*sin(x[4])*sin(x[3])+(3*mp*r*cos(x[7])*x[15]+m*x[9])*sin(x[4])*cos(x[3]))*x[11]+((3*mp*r*cos(x[7])*x[15]+m*x[9])*cos(x[4])*x[12]+(mp*r*cos(x[7])*k[15]-mp*r*sin(x[7])*pow(x[15],2))*sin(x[4]))*sin(x[3])+((3*mp*r*cos(x[7])*sin(x[6])*x[14]+3*mp*r*sin(x[7])*x[15]*cos(x[6])+m*x[10])*cos(x[4])*x[12]+(mp*r*cos(x[7])*sin(x[6])*k[14]+mp*r*cos(x[7])*cos(x[6])*pow(x[14],2)-2*mp*r*sin(x[7])*x[15]*sin(x[6])*x[14]+(mp*r*sin(x[7])*k[15]+mp*r*cos(x[7])*pow(x[15],2))*cos(x[6]))*sin(x[4]))*cos(x[3])+((-3*mp*r*cos(x[7])*cos(x[6])*x[14])+3*mp*r*sin(x[7])*x[15]*sin(x[6])-m*x[8])*sin(x[4])*x[12]+(mp*r*cos(x[7])*cos(x[6])*k[14]-mp*r*cos(x[7])*sin(x[6])*pow(x[14],2)-2*mp*r*sin(x[7])*x[15]*cos(x[6])*x[14]+((-mp*r*sin(x[7])*k[15])-mp*r*cos(x[7])*pow(x[15],2))*sin(x[6]))*cos(x[4]))*sin(x[5])+((((-3*mp*r*cos(x[7])*x[15])-m*x[9])*sin(x[3])+((-3*mp*r*cos(x[7])*sin(x[6])*x[14])-3*mp*r*sin(x[7])*x[15]*cos(x[6])-m*x[10])*cos(x[3]))*x[11]+((-mp*r*cos(x[7])*sin(x[6])*k[14])-mp*r*cos(x[7])*cos(x[6])*pow(x[14],2)+2*mp*r*sin(x[7])*x[15]*sin(x[6])*x[14]+((-mp*r*sin(x[7])*k[15])-mp*r*cos(x[7])*pow(x[15],2))*cos(x[6]))*sin(x[3])+(mp*r*cos(x[7])*k[15]-mp*r*sin(x[7])*pow(x[15],2))*cos(x[3]))*cos(x[5])-y_dir*U1);
	k[10] = pow(mp+m,-1)*(((3*mp*r*cos(x[7])*sin(x[6])*x[14]+3*mp*r*sin(x[7])*x[15]*cos(x[6])+m*x[10])*cos(x[4])*sin(x[3])+((-3*mp*r*cos(x[7])*x[15])-m*x[9])*cos(x[4])*cos(x[3]))*x[11]+((3*mp*r*cos(x[7])*x[15]+m*x[9])*sin(x[4])*x[12]+(mp*r*sin(x[7])*pow(x[15],2)-mp*r*cos(x[7])*k[15])*cos(x[4]))*sin(x[3])+((3*mp*r*cos(x[7])*sin(x[6])*x[14]+3*mp*r*sin(x[7])*x[15]*cos(x[6])+m*x[10])*sin(x[4])*x[12]+((-mp*r*cos(x[7])*sin(x[6])*k[14])-mp*r*cos(x[7])*cos(x[6])*pow(x[14],2)+2*mp*r*sin(x[7])*x[15]*sin(x[6])*x[14]+((-mp*r*sin(x[7])*k[15])-mp*r*cos(x[7])*pow(x[15],2))*cos(x[6]))*cos(x[4]))*cos(x[3])+(3*mp*r*cos(x[7])*cos(x[6])*x[14]-3*mp*r*sin(x[7])*x[15]*sin(x[6])+m*x[8])*cos(x[4])*x[12]+(mp*r*cos(x[7])*cos(x[6])*k[14]-mp*r*cos(x[7])*sin(x[6])*pow(x[14],2)-2*mp*r*sin(x[7])*x[15]*cos(x[6])*x[14]+((-mp*r*sin(x[7])*k[15])-mp*r*cos(x[7])*pow(x[15],2))*sin(x[6]))*sin(x[4])-g*mp-g*m+z_dir*U1);

	k[11] = -pow(Ixx,-1)*pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-1)*(((Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*cos(x[4])*pow(sin(x[4]),2)+(Iyy*pow(Izz,2)-pow(Iyy,2)*Izz)*pow(cos(x[4]),3))*cos(x[3])*sin(x[3])*pow(x[11],2)+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*cos(x[4])*sin(x[4])*cos(x[3])*sin(x[3])*x[9]+((Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*pow(sin(x[4]),2)+(2*Iyy*pow(Izz,2)-2*pow(Iyy,2)*Izz)*pow(cos(x[4]),2))*x[10]*pow(cos(x[3]),2)+(((pow(Ixx,2)-Ixx*Iyy)*Izz-Ixx*pow(Izz,2))*pow(sin(x[4]),2)+((pow(Iyy,2)-Ixx*Iyy)*Izz-Iyy*pow(Izz,2))*pow(cos(x[4]),2))*x[10])*x[11]+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*sin(x[4])*x[10]*pow(cos(x[3]),2)+(Ixx*pow(Izz,2)+((-Ixx*Iyy)-pow(Ixx,2))*Izz)*sin(x[4])*x[10])*x[9]+((pow(Iyy,2)*Izz-Iyy*pow(Izz,2))*cos(x[4])*pow(x[10],2)*cos(x[3])-Ixx*Izz*U3*sin(x[4]))*sin(x[3])-Ixx*Iyy*U4*sin(x[4])*cos(x[3])-Iyy*Izz*U2*cos(x[4]));
	k[12] = pow(Iyy,-1)*pow(Izz,-1)*(((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*pow(sin(x[3]),2)+(Ixx*Izz-pow(Izz,2))*cos(x[4])*sin(x[4]))*pow(x[11],2)+((((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*pow(sin(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*cos(x[4]))*x[9]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*sin(x[4])*x[10]*cos(x[3])*sin(x[3]))*x[11]+((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*x[10]*cos(x[3])*sin(x[3])*x[9]-Iyy*U4*sin(x[3])+Izz*U3*cos(x[3]));
	k[13] = -pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-1)*((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*cos(x[3])*sin(x[3])*pow(x[11],2)+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*cos(x[3])*sin(x[3])*x[9]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*sin(x[4])*x[10]*pow(cos(x[3]),2)+((Ixx-Iyy)*Izz-pow(Izz,2))*sin(x[4])*x[10])*x[11]+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*x[10]*pow(cos(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*x[10])*x[9]-Izz*U3*sin(x[3])-Iyy*U4*cos(x[3]));

	k[14] = 0.0;
	k[15] = 0.0;
}

int main(void)
{
	int i, j;

	double t = 0.0;
	double x[X_NUM] = {X_POS, Y_POS, Z_POS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};	
	double c[RK4_SIZE] = {0.0, 0.5, 0.5, 1.0};
	double **state;
	double **k;

	state = malloc(sizeof(double *) * RK4_SIZE);
	k = malloc(sizeof(double *) * RK4_SIZE);
	i = 0;
	while(i < RK4_SIZE)
	{
		state[i] = malloc(sizeof(double) * X_NUM);
		k[i] = malloc(sizeof(double) * X_NUM);
		j = 0;
		while (j < X_NUM)
		{
			state[i][j] = 0.0;
			k[i][j] = 0.0;
			j++;
		}
		i++;
	}

	while (t <= 10.0)
	{
		printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15]);
		
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
			func(t, state[i], k[i]);
			i++;
		}
		j = 0;
		while (j < X_NUM)
		{
			x[j] += ((k[0][j] + 2.0*k[1][j] + 2.0*k[2][j] + k[3][j])*(dt))/6.0;
			j++;
		}
		t += dt;
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
