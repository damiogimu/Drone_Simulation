#include "sim.h"

#define Z_DES 2.0

void any_traj(double t, t_desire *des)
{
	if (t <= Z_RISE_T)
	{
		des->xd = 0.0;
		des->dxd = 0.0;
		des->ddxd = 0.0;
		des->yd = 0.0;
		des->dyd = 0.0;
		des->ddyd = 0.0;
	}
	else
	{
		t -= Z_RISE_T;
		des->xd = 0.1*t*t;
		des->dxd = 0.2*t;
		des->ddxd = 0.2;
		des->yd = 0.1*t*t;
		des->dyd = 0.2*t;
		des->ddyd = 0.2;
	}
	des->zd = Z_DES;
	des->dzd = 0.0;
	des->ddzd = 0.0;
	des->zd += 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}

void noncontrol_traj_vt(double t, t_desire *des)
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
	des->zd -= 1.0e-8;
	des->cd = 0.0;
	des->dcd = 0.0;
	des->ddcd = 0.0;
}

void noncontrol_traj_xt(double t, t_desire *des)
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
