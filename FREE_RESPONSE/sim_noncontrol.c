#include "my_header.h"

#define SCALE 1.0

void func(double t, double *x, double *k)
{
	double U1, U2, U3, U4;
	double x_dir = (cos(x[3])*sin(x[4])*cos(x[5])) + (sin(x[3])*sin(x[5]));
	double y_dir = (cos(x[3])*sin(x[4])*sin(x[5])) - (sin(x[3])*cos(x[5]));
	double z_dir = cos(x[3])*cos(x[4]);
	
	U1 = (m+mp)*g;
	U2 = 0.001;
	U3 = 0.001;
	U4 = 0.0;

	k[0] = x[8];
	k[1] = x[9];
	k[2] = x[10];
	k[3] = x[11];
	k[4] = x[12];
	k[5] = x[13];
	k[6] = x[14];
	k[7] = x[15];

	k[8] = pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),3)*sin(x[6])*pow(x[14],2)-x_dir*U1*mp*pow(cos(x[7]),2)*pow(sin(x[6]),2)+(z_dir*U1*mp*pow(cos(x[7]),2)*cos(x[6])+m*mp*r*cos(x[7])*pow(x[15],2)-y_dir*U1*mp*cos(x[7])*sin(x[7]))*sin(x[6])+x_dir*U1*mp+x_dir*U1*m);
	k[9] = pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),2)*sin(x[7])*pow(x[14],2)-x_dir*U1*mp*cos(x[7])*sin(x[7])*sin(x[6])+z_dir*U1*mp*cos(x[7])*sin(x[7])*cos(x[6])+m*mp*r*sin(x[7])*pow(x[15],2)+y_dir*U1*mp*pow(cos(x[7]),2)+y_dir*U1*m);
	k[10] = -pow(m*mp+pow(m,2),-1)*(m*mp*r*pow(cos(x[7]),3)*cos(x[6])*pow(x[14],2)-x_dir*U1*mp*pow(cos(x[7]),2)*cos(x[6])*sin(x[6])+z_dir*U1*mp*pow(cos(x[7]),2)*pow(cos(x[6]),2)+(m*mp*r*cos(x[7])*pow(x[15],2)-y_dir*U1*mp*cos(x[7])*sin(x[7]))*cos(x[6])+(g*m-z_dir*U1)*mp+g*pow(m,2)-z_dir*U1*m);
/*
	k[11] = -pow(Ixx,-1)*pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-2)*((((Iyy-Ixx)*pow(Izz,2)+(pow(Ixx,2)-pow(Iyy,2))*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),4)+(Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*pow(cos(x[4]),2))*cos(x[3])*sin(x[3])*pow(x[13],2)+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),2)*sin(x[4])*cos(x[3])*sin(x[3])*x[11]+(((2*Iyy-Ixx)*pow(Izz,2)+(pow(Ixx,2)-2*pow(Iyy,2))*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*pow(cos(x[4]),3)+(Ixx*pow(Izz,2)-pow(Ixx,2)*Izz-Ixx*pow(Iyy,2)+pow(Ixx,2)*Iyy)*cos(x[4]))*x[12]*pow(cos(x[3]),2)+(((Ixx-Iyy)*pow(Izz,2)+(pow(Iyy,2)-pow(Ixx,2))*Izz)*pow(cos(x[4]),3)+((pow(Ixx,2)-Ixx*Iyy)*Izz-Ixx*pow(Izz,2))*cos(x[4]))*x[12])*x[13]+(((-Ixx*pow(Izz,2))+pow(Ixx,2)*Izz+Ixx*pow(Iyy,2)-pow(Ixx,2)*Iyy)*cos(x[4])*sin(x[4])*x[12]*pow(cos(x[3]),2)+(Ixx*pow(Izz,2)+((-Ixx*Iyy)-pow(Ixx,2))*Izz)*cos(x[4])*sin(x[4])*x[12])*x[11]+((pow(Iyy,2)*Izz-Iyy*pow(Izz,2))*pow(cos(x[4]),2)*pow(x[12],2)+(Ixx*Iyy-Ixx*Izz)*U3*cos(x[4])*sin(x[4]))*cos(x[3])*sin(x[3])+((Ixx*Izz-Ixx*Iyy)*U4*sin(x[4])+(Ixx*Iyy-Ixx*Izz)*U2*pow(cos(x[4]),2)+(Ixx*Izz-Ixx*Iyy)*U2)*pow(cos(x[3]),2)-Ixx*Izz*U4*sin(x[4])+(Ixx-Iyy)*Izz*U2*pow(cos(x[4]),2)-Ixx*Izz*U2);
	k[12] = pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-1)*(((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*pow(cos(x[4]),2)*sin(x[4])*pow(sin(x[3]),2)+(Ixx*Izz-pow(Izz,2))*pow(cos(x[4]),2)*sin(x[4]))*pow(x[13],2)+((((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*pow(cos(x[4]),2)*pow(sin(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*pow(cos(x[4]),2))*x[11]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*x[12]*cos(x[3])*sin(x[3]))*x[13]+((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*x[12]*cos(x[3])*sin(x[3])*x[11]+(Iyy-Izz)*U3*cos(x[4])*pow(sin(x[3]),2)+((Izz-Iyy)*U2*sin(x[4])+(Izz-Iyy)*U4)*cos(x[3])*sin(x[3])+Izz*U3*cos(x[4]));
	k[13] = -pow(Iyy,-1)*pow(Izz,-1)*pow(cos(x[4]),-2)*((pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*pow(cos(x[4]),2)*sin(x[4])*cos(x[3])*sin(x[3])*pow(x[13],2)+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*pow(cos(x[4]),2)*cos(x[3])*sin(x[3])*x[11]+(pow(Izz,2)-Ixx*Izz-pow(Iyy,2)+Ixx*Iyy)*cos(x[4])*sin(x[4])*x[12]*pow(cos(x[3]),2)+((Ixx-Iyy)*Izz-pow(Izz,2))*cos(x[4])*sin(x[4])*x[12])*x[13]+(((-pow(Izz,2))+Ixx*Izz+pow(Iyy,2)-Ixx*Iyy)*cos(x[4])*x[12]*pow(cos(x[3]),2)+(pow(Izz,2)+((-Iyy)-Ixx)*Izz)*cos(x[4])*x[12])*x[11]+(Iyy-Izz)*U3*cos(x[4])*cos(x[3])*sin(x[3])+((Izz-Iyy)*U2*sin(x[4])+(Izz-Iyy)*U4)*pow(cos(x[3]),2)-Izz*U2*sin(x[4])-Izz*U4);
*/
	k[11] = (x[12]*x[13]*(Iyy - Izz))/Ixx + U2/Ixx;
	k[12] = (x[11]*x[13]*(Izz - Ixx))/Ixx + U3/Iyy;
	k[13] = (x[11]*x[12]*(Ixx - Iyy))/Ixx + U4/Izz;
	
	k[14] = pow(m,-1)*pow(r,-1)*pow(cos(x[7]),-1)*(2*m*r*sin(x[7])*x[15]*x[14]-z_dir*U1*sin(x[6])-x_dir*U1*cos(x[6]));
	k[15] = -pow(m,-1)*pow(r,-1)*(m*r*cos(x[7])*sin(x[7])*pow(x[14],2)-x_dir*U1*sin(x[7])*sin(x[6])+z_dir*U1*sin(x[7])*cos(x[6])+y_dir*U1*cos(x[7]));
}

int main(void)
{
	FILE *fd[FD_NUM];
	char files[FD_NUM][100]={DATA_FILE1, DATA_FILE2, DATA_FILE3, DATA_FILE4, DATA_FILE5};
	int i, j;
	double t = 0.0;
	double x[X_NUM] = {INIT_X, INIT_Y, INIT_Z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double c[RK4_SIZE] = {0.0, 0.5, 0.5, 1.0};
	double **state;
	double **k;
	double ddx, ddy, ddz;
	double x_p, y_p, z_p;
	t_rotor rot;

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
	while (i < FD_NUM)
	{
		fd[i] = fopen(files[i], "w");
		i++;
	}

	while (t <= TIME)
	{
//----------  output data  -----------//
		x_p = x[0] + r*cos(x[7])*sin(x[6]);
		y_p = x[1] + r*sin(x[7]);
		z_p = x[2] - r*cos(x[7])*cos(x[6]);
		fprintf(fd[0], "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ", t, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15], x_p, y_p, z_p);
		if ((int)((t*1000000.0)+0.5)%50000==0)
		{
			fprintf(fd[1], "%f %f %f %f\n", t,x[0],x[1],x[2]);
			fprintf(fd[1], "%f %f %f %f\n", t,x_p,y_p,z_p);
			fprintf(fd[1], "\n\n");
			rot.l_x1 = x[0]+(SCALE*l)*cos(x[5])*cos(x[5]);
			rot.l_y1 = x[1]+(SCALE*l)*sin(x[5])*cos(x[5]);
			rot.l_z1 = x[2] -(SCALE*l)*tan(x[4]);
			fprintf(fd[3], "%f %f %f %f\n", t,rot.l_x1,rot.l_y1,rot.l_z1);
			rot.l_x2 = x[0]-(SCALE*l)*cos(x[5])*cos(x[5]);
			rot.l_y2 = x[1]-(SCALE*l*cos(x[5])*sin(x[5]));
			rot.l_z2 = x[2] + (SCALE*l)*tan(x[4]);
			fprintf(fd[3], "%f %f %f %f\n", t,rot.l_x2,rot.l_y2,rot.l_z2);
			fprintf(fd[3], "\n\n");
			rot.l_x1 = x[0]+(SCALE*l*cos(x[5])*sin(x[5]));
			rot.l_y1 = x[1]-(SCALE*l)*cos(x[5])*cos(x[5]);
			rot.l_z1 = x[2] - (SCALE*l)*tan(x[3]);
			fprintf(fd[4], "%f %f %f %f\n", t,rot.l_x1,rot.l_y1,rot.l_z1);
			rot.l_x2 = x[0]-(SCALE*l*sin(x[5])*cos(x[5]));
			rot.l_y2 = x[1]+(SCALE*l)*cos(x[5])*cos(x[5]);
			rot.l_z2 = x[2] + (SCALE*l)*tan(x[3]);
			fprintf(fd[4], "%f %f %f %f\n", t,rot.l_x2,rot.l_y2,rot.l_z2);
			fprintf(fd[4], "\n\n");
			fprintf(fd[2], "%f %f %f %f %f %f %f %f %f %f\n", t,x[0],x[1],x[2],x[3],x[4],x[5],x_p,y_p,z_p);
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
			func(t, state[i], k[i]);
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
		fprintf(fd[0], "%f %f %f\n", ddx, ddy, ddz);
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
