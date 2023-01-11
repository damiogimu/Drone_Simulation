#include "my_header.h"

#define DDX 0.1
#define DDY 0.1
#define DDZ 0.0

void func(double t, double *x, double *k, double THRUST)
{
	double U1, U2, U3, U4;
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
	U1 = m*g;
	U2 = 0.0;
	U3 = 0.0;
	U4 = 0.0;

	k[0] = x[6];
	k[1] = x[7];
	k[2] = x[8];
	k[3] = x[9];
	k[4] = x[10];
	k[5] = x[11];

	k[6] = pow(m,-1)*(U1*sin(x[3])*sin(x[5])+U1*sin(x[4])*cos(x[3])*cos(x[5]));
	k[7] = pow(m,-1)*(U1*sin(x[4])*cos(x[3])*sin(x[5])-U1*sin(x[3])*cos(x[5]));
	k[8] = pow(m,-1)*(U1*cos(x[4])*cos(x[3])-g*m);
	k[9] = -pow(Ixx,-1)*pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-2)*((((Iyy-Ixx)*pow(Izz,2)+(pow(Ixx,2)-pow(Iyy,2))*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),4)+(Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*pow(cos(x[4]),2))*cos(x[3])*sin(x[3])*pow(x[11],2)+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),2)*sin(x[4])*cos(x[3])*sin(x[3])*x[9]+(((2*Iyy-Ixx)*pow(Izz,2)+(pow(Ixx,2)-2*pow(Iyy,2))*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),3)+(Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*cos(x[4]))*x[10]*pow(cos(x[3]),2)+(((Ixx-Iyy)*pow(Izz,2)+(pow(Iyy,2)-pow(Ixx,2))*Izz)*pow(cos(x[4]),3)+((pow(Ixx,2)-Ixx*Iyy)*Izz-Ixx*pow(Izz,2))*cos(x[4]))*x[10])*x[11]+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*cos(x[4])*sin(x[4])*x[10]*pow(cos(x[3]),2)+(Ixx*pow(Izz,2)+((-Ixx*Iyy)-pow(Ixx,2))*Izz)*cos(x[4])*sin(x[4])*x[10])*x[9]+((pow(Iyy,2)*Izz-Iyy*pow(Izz,2))*pow(cos(x[4]),2)*pow(x[10],2)+(Ixx*Iyy-Ixx*Izz)*U3*cos(x[4])*sin(x[4]))*cos(x[3])*sin(x[3])+((Ixx*Izz-Ixx*Iyy)*U4*sin(x[4])+(Ixx*Iyy-Ixx*Izz)*U2*pow(cos(x[4]),2)+(Ixx*Izz-Ixx*Iyy)*U2)*pow(cos(x[3]),2)-Ixx*Izz*U4*sin(x[4])+(Ixx-Iyy)*Izz*U2*pow(cos(x[4]),2)-Ixx*Izz*U2);
	k[10] = pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-1)*(((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*pow(cos(x[4]),2)*sin(x[4])*pow(sin(x[3]),2)+(Ixx*Izz-pow(Izz,2))*pow(cos(x[4]),2)*sin(x[4]))*pow(x[11],2)+((((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*pow(cos(x[4]),2)*pow(sin(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*pow(cos(x[4]),2))*x[9]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*x[10]*cos(x[3])*sin(x[3]))*x[11]+((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*x[10]*cos(x[3])*sin(x[3])*x[9]+(Iyy-Izz)*U3*cos(x[4])*pow(sin(x[3]),2)+((Izz-Iyy)*U2*sin(x[4])+(Izz-Iyy)*U4)*cos(x[3])*sin(x[3])+Izz*U3*cos(x[4]));
	k[11] = -pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-2)*((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*pow(cos(x[4]),2)*sin(x[4])*cos(x[3])*sin(x[3])*pow(x[11],2)+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*pow(cos(x[4]),2)*cos(x[3])*sin(x[3])*x[9]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*x[10]*pow(cos(x[3]),2)+((Ixx-Iyy)*Izz-pow(Izz,2))*cos(x[4])*sin(x[4])*x[10])*x[11]+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*x[10]*pow(cos(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*cos(x[4])*x[10])*x[9]+(Iyy-Izz)*U3*cos(x[4])*cos(x[3])*sin(x[3])+((Izz-Iyy)*U2*sin(x[4])+(Izz-Iyy)*U4)*pow(cos(x[3]),2)-Izz*U2*sin(x[4])-Izz*U4);
}

int main(void)
{
	int i, j;

/*
	double THETA = atan(DDX/(DDZ+g));
	double PHI = -atan(DDY/(DDX*sin(THETA) + DDZ*cos(THETA) + g*cos(THETA)));
	double THRUST = m*(DDZ+g)/(cos(THETA)*cos(PHI));
*/
	double t = 0.0;
	double x[X_NUM] = {X_POS, Y_POS, Z_POS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};	
	double c[RK4_SIZE] = {0.0, 0.5, 0.5, 1.0};
	double **state;
	double **k;
	double ddx, ddy, ddz;

	state = malloc(sizeof(double *) * RK4_SIZE);
	k = malloc(sizeof(double *) * RK4_SIZE);
	i = 0;
	while(i < RK4_SIZE)
	{
		state[i] = malloc(sizeof(double) * X_NUM);
		k[i] = malloc(sizeof(double) * X_NUM);
		i++;
	}

	while (t <= TIME)
	{
		printf("%f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]);
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
			func(t, state[i], k[i], 0.0);
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
