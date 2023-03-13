#include "sim.h"

#define Z_DES 2.0

void any_traj(double t, t_desire *des)
{
	des->xd = 0.0;
	des->dxd = 0.0;
	des->ddxd = 0.0;
	des->yd = 0.0;
	des->dyd = 0.0;
	des->ddyd = 0.0;
	des->zd = 2.0;
	des->dzd = Z_DES;
	des->ddzd = 0.0;
	des->zd += 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}

void noncontrol_traj1(double t, t_desire *des)
{
	double T;
	T = vmax/amax;
	des->ACC_T = T;
	
	if (t <= Z_RISE_T)
	{
		des->xd = 0.0;
		des->dxd = 0.0;
		des->ddxd = 0.0;
		des->yd = 0.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	else if (t <= Z_RISE_T+T)
	{
		t -= Z_RISE_T;
		des->xd = (amax*t*t)/2.0;
		des->dxd = amax*t;
		des->ddxd = amax;
		des->yd = (amax*t*t)/2.0;
		des->dyd = amax*t;
		des->ddyd = amax;
	}
	else
	{
		t -= Z_RISE_T+T;
		des->xd = amax*T*t + (amax*T*T)/2.0;
		des->dxd = amax*T;
		des->ddxd = 0.0;
		des->yd = amax*T*t + (amax*T*T)/2.0;
		des->dyd = amax*T;
		des->ddyd = 0.0;
	}
	des->zd = Z_DES;
	des->dzd = 0.0;
	des->ddzd = 0.0;
	des->zd += 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}

void controled_traj1(double t, double f, t_desire *des)
{
	double a, T, T1, dT;
	
	T = vmax/amax;
	T1 = T/2.0;
	a = ceil((f*T-1.0)/2.0);
	dT = (a+0.5)/f - T/2.0;
	des->ACC_T = T+dT;
	
	if (t <= Z_RISE_T)
	{
		des->xd = 0.0;
		des->dxd = 0.0;
		des->ddxd = 0.0;
		des->yd = 0.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	else if (t < Z_RISE_T+T1)
	{
		t -= Z_RISE_T;
		des->xd = (amax*t*t)/2.0;
		des->dxd = amax*t;
		des->ddxd = amax;
		des->yd = (amax*t*t)/2.0;
		des->dyd = amax*t;
		des->ddyd = amax;
	}
	else if (t < Z_RISE_T+T1+dT)
	{
		t -= Z_RISE_T+T1;
		des->xd = amax*T1*t + (amax*T1*T1)/2.0;
		des->dxd = amax*T1;
		des->ddxd = 0.0;
		des->yd = amax*T1*t + (amax*T1*T1)/2.0;
		des->dyd = amax*T1;
		des->ddyd = 0.0;
	}
	else if (t < Z_RISE_T+T+dT)
	{
		t -= Z_RISE_T+T1+dT;
		des->xd = (amax*t*t)/2.0 + amax*T1*t + amax*T1*dT + (amax*T1*T1)/2.0;
		des->dxd = amax*t + amax*T1;
		des->ddxd = amax;
		des->yd = (amax*t*t)/2.0 + amax*T1*t + amax*T1*dT + (amax*T1*T1)/2.0;
		des->dyd = amax*t + amax*T1;
		des->ddyd = amax;
	}
	else
	{
		t -= Z_RISE_T+T+dT;
		des->xd = 2.0*amax*T1*t + (amax*T1*T1)/2.0 + amax*T1*T1 + amax*T1*dT + (amax*T1*T1)/2.0;
		des->dxd = 2.0*amax*T1;
		des->ddxd = 0.0;
		des->yd = 2.0*amax*T1*t + (amax*T1*T1)/2.0 + amax*T1*T1 + amax*T1*dT + (amax*T1*T1)/2.0;
		des->dyd = 2.0*amax*T1;
		des->ddyd = 0.0;
	}
	des->zd = Z_DES;
	des->dzd = 0.0;
	des->ddzd = 0.0;
	des->zd += 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}

void noncontrol_traj2(double t, t_desire *des)
{
	double tx, Tx, T2x;
	double ty, Ty, T2y;

	tx = t;
	ty = t;
	Tx = sqrt(Xt/amax);
	Ty = sqrt(Yt/amax);
	T2x = 2.0*Tx;
	T2y = 2.0*Ty;
	des->ACC_T = T2x;

	if (tx <= Z_RISE_T)
	{
		des->xd = 0.0;
		des->dxd = 0.0;
		des->ddxd = 0.0;
	}
	else if (tx <= Z_RISE_T+Tx)
	{
		tx -= Z_RISE_T;
		des->xd = (amax*tx*tx)/2.0;
		des->dxd = amax*tx;
		des->ddxd = amax;
	}
	else if (tx <= Z_RISE_T+T2x)
	{
		tx -= Z_RISE_T+Tx;
		des->xd = -(amax*tx*tx)/2.0 + amax*Tx*tx + (amax*Tx*Tx)/2.0;
		des->dxd = -amax*tx + amax*Tx;
		des->ddxd = -amax;
	}
	else
	{
		tx -= Z_RISE_T+T2x;
		des->xd = -(amax*Tx*Tx)/2.0 + amax*Tx*Tx + (amax*Tx*Tx)/2.0;
		des->dxd = 0.0;
		des->ddxd = 0.0;
	}

	if (ty <= Z_RISE_T)
	{
		des->yd = 0.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	else if (ty <= Z_RISE_T+Ty)
	{
		ty -= Z_RISE_T;
		des->yd = (amax*ty*ty)/2.0;
		des->dyd = amax*ty;
		des->ddyd = amax;
	}
	else if (ty <= Z_RISE_T+T2y)
	{
		ty -= Z_RISE_T+Ty;
		des->yd = -(amax*ty*ty)/2.0 + amax*Ty*ty + (amax*Ty*Ty)/2.0;
		des->dyd = -amax*ty + amax*Ty;
		des->ddyd = -amax;
	}
	else
	{
		ty -= Z_RISE_T+T2y;
		des->yd = -(amax*Ty*Ty)/2.0 + amax*Ty*Ty + (amax*Ty*Ty)/2.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	des->zd = Z_DES;
	des->dzd = 0.0;
	des->ddzd = 0.0;
	des->zd += 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}

void controled_traj2(double t, double f, t_desire *des)
{
	double a, tx, Tx, T1x, dTx;
	double b, ty, Ty, T1y, dTy;
	
	tx = ty = t;
	a = ceil(f*(sqrt(Xt/amax)));
	if (a == 0.0)
		a = 1.0;
	Tx = (2*f*Xt)/(a*amax);
	T1x = Tx/2.0;
	dTx = a/f - Tx/2.0;
	b = ceil(f*(sqrt(Yt/amax)));
	if (b == 0.0)
		b = 1.0;
	Ty = (2*f*Yt)/(b*amax);
	T1y = Ty/2.0;
	dTy = b/f - Ty/2.0;
	des->ACC_T = Tx+dTx;

	if (tx <= Z_RISE_T)
	{
		des->xd = 0.0;
		des->dxd = 0.0;
		des->ddxd = 0.0;
	}
	else if (tx <= Z_RISE_T+T1x)
	{
		tx -= Z_RISE_T;
		des->xd = (amax*tx*tx)/2.0;
		des->dxd = amax*tx;
		des->ddxd = amax;
	}
	else if (tx <= Z_RISE_T+T1x+dTx)
	{
		tx -= Z_RISE_T+T1x;
		des->xd = amax*T1x*tx + (amax*T1x*T1x)/2.0;
		des->dxd = amax*T1x;
		des->ddxd = 0.0;
	}
	else if (tx <= Z_RISE_T+Tx+dTx)
	{
		tx -= Z_RISE_T+T1x+dTx;
		des->xd = (-amax*tx*tx)/2.0 + amax*T1x*tx + amax*T1x*dTx + (amax*T1x*T1x)/2.0;
		des->dxd = -amax*tx + amax*T1x;
		des->ddxd = -amax;
	}
	else
	{
		des->xd = (-amax*T1x*T1x)/2.0 + amax*T1x*T1x + amax*T1x*dTx + (amax*T1x*T1x)/2.0;
		des->dxd = -amax*T1x + amax*T1x;
		des->ddxd = 0.0;
	}

	if (ty <= Z_RISE_T)
	{
		des->yd = 0.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	else if (ty <= Z_RISE_T+T1y)
	{
		ty -= Z_RISE_T;
		des->yd = (amax*ty*ty)/2.0;
		des->dyd = amax*ty;
		des->ddyd = amax;
	}
	else if (ty <= Z_RISE_T+T1y+dTy)
	{
		ty -= Z_RISE_T+T1y;
		des->yd = amax*T1y*ty + (amax*T1y*T1y)/2.0;
		des->dyd = amax*T1y;
		des->ddyd = 0.0;
	}
	else if (ty <= Z_RISE_T+Ty+dTy)
	{
		ty -= Z_RISE_T+T1y+dTy;
		des->yd = (-amax*ty*ty)/2.0 + amax*T1y*ty + amax*T1y*dTy + (amax*T1y*T1y)/2.0;
		des->dyd = -amax*ty + amax*T1y;
		des->ddyd = -amax;
	}
	else
	{
		des->yd = (-amax*T1y*T1y)/2.0 + amax*T1y*T1y + amax*T1y*dTy + (amax*T1y*T1y)/2.0;
		des->dyd = -amax*T1y + amax*T1y;
		des->ddyd = 0.0;
	}
	des->zd = Z_DES;
	des->dzd = 0.0;
	des->ddzd = 0.0;
	des->zd += 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}
