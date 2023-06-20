#include"my_header.h"

void set_desire_path1(double t, double f, t_desire *des)
{
	double a, T, T1, dT;

	T = vmax/amax;
	T1 = T/2.0;
	a = ceil((f*T-1.0)/2.0);
	dT = (a+0.5)/f - T/2.0;

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
		des->yd = 0.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	else if (t < Z_RISE_T+T1+dT)
	{
		t -= Z_RISE_T+T1;
		des->xd = amax*T1*t + (amax*T1*T1)/2.0;
		des->dxd = amax*T1;
		des->ddxd = 0.0;
		des->yd = 0.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	else if (t < Z_RISE_T+T+dT)
	{
		t -= Z_RISE_T+T1+dT;
		des->xd = (amax*t*t)/2.0 + amax*T1*t + amax*T1*dT + (amax*T1*T1)/2.0;
		des->dxd = amax*t + amax*T1;
		des->ddxd = amax;
		des->yd = 0.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	else
	{
		t -= Z_RISE_T+T+dT;
		des->xd = 2.0*amax*T1*t + (amax*T1*T1)/2.0 + amax*T1*T1 + amax*T1*dT + (amax*T1*T1)/2.0;
		des->dxd = 2.0*amax*T1;
		des->ddxd = 0.0;
		des->yd = 0.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	des->zd = 2.0;
	des->dzd = 0.0;
	des->ddzd = 0.0;
	des->zd += 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}

void set_desire_path2(double t, double f, t_desire *des)
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

	des->zd = 2.0;
	des->dzd = 0.0;
	des->ddzd = 0.0;
	des->zd += 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}

void vertical_traj(double t, t_desire *des)
{
	des->xd = 0.0;
	des->dxd = 0.0;
	des->ddxd = 0.0;
	des->yd = 0.0;
	des->dyd = 0.0;
	des->ddyd = 0.0;
	des->zd = 2.0;
	des->dzd = 0.0;
	des->ddzd = 0.0;
	des->zd += 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}
