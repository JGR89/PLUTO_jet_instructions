#include "pluto.h"

#define min_tracer_time  0.0 // time of first tracer injection
#define max_tracer_time  4.0 // time of final tracer injection
#define inject_frac  0.05 // fraction of time tracer is

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
    int t;
    double r, cs, Temp;
    r = x1;
    
   // v[RHO] = pow((1+(x1/g_inputParam[R_CORE])*(x1/g_inputParam[R_CORE])),-1.5*g_inputParam[B_EXPONENT]); /* King profile of Krause 2005 */
    cs = 1;
    v[RHO] = pow(1 + pow(r/g_inputParam[R_CORE],2),-1.5*g_inputParam[B_EXPONENT]);
    g_gamma = 5.0/3.0;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[PRS]  = (cs/g_gamma)*v[RHO];
    if (NTRACER > 0)
    {
        v[TRC + 0] = 0.0;   /* jet tracer */
        
        for (t = 1; t < NTRACER; t++) //so start at t, and incriment t until t = NTRACER
        {
            //printf("tracer1 %d\n", t);
            v[TRC + t] = 0.0;   /* tracer i */
        }
    }
}
/* **************************************************************** */

void Analysis (const Data *d, Grid *grid)

/* **************************************************************** */
{
}
/* **************************************************************** */

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)

/* ***********************************************************************
 *  Assign user-defined boundary conditions in the lower boundary ghost  *
 *  zones.  Top hat conical using L1abc, see Krause+2012.                *
 *                                                                       *
 *************************************************************************/
{
    int   i, j, k, nv, t;
    double  *x1, *x2, *x3;
    double cs, Tj, vj[NVAR];
    double  *r, *theta, L1b, Omega, theta_rad;
    
    theta_rad = g_inputParam[H_OPEN_ANG] * CONST_PI / 180.0;
    Omega= 2.0*CONST_PI*(1.0-cos(theta_rad));
    L1b  = 0.5/sqrt(Omega);
    r = grid[IDIR].xgc;  /* -- array pointer to x1 coordinate -- */
    theta = grid[JDIR].xgc;  /* -- array pointer to x2 coordinate -- */
    /* x3 = grid[KDIR].xgc;   -- array pointer to x3 coordinate -- */
    //printf("L1b is %e\n theta_rad is %e\n Omega is %e\n input angle is %e\n ",L1b, theta_rad, Omega, g_inputParam[H_OPEN_ANG]);
    vj[RHO] = pow(L1b,2);
    //printf("Calculated rhojet is %e. I input %e\n",pow(L1b,2),vj[RHO]);
    vj[PRS] = 1.0/g_gamma;         /* -- Pressure-matched jet -- */
    vj[VX1] = g_inputParam[MACH_EXT];  /* -- Sound speed is one in this case -- */
    
    if (side == 0)
    {
        TOT_LOOP(k,j,i)
        {
            if (d->Vc[PRS][k][j][i] < 1.e-6) /* set min pressure */
            {
                d->Vc[PRS][k][j][i] = 1.e-6;
            }
        }
    }
    
    if (side == X1_BEG)
    {     /* -- select the boundary side -- */
        BOX_LOOP(box,k,j,i)
        {   /* -- Loop over boundary zones -- */
            if (theta[j] <= theta_rad) /* -- set jet values for r <= theta_rad -- */
            {
                d->Vc[RHO][k][j][i] = vj[RHO] / pow(r[i],2);
                d->Vc[VX2][k][j][i] = 0.0;
                d->Vc[VX1][k][j][i] = vj[VX1];
                d->Vc[PRS][k][j][i] = vj[PRS];
                if (NTRACER > 0)
                {
                    d->Vc[TRC][k][j][i] = 1.0; /* tracer particle over jet */
                    
                    for (t = 1; t < NTRACER; t++) /* -- inject tracer particles --*/
                    {
                        if (min_tracer_time + (t - 1)*(max_tracer_time - min_tracer_time)/((double) NTRACER - 1) <= g_time && g_time < min_tracer_time + (t - 1 + inject_frac)*(max_tracer_time - min_tracer_time)/((double) NTRACER - 1))
                        {
                            d->Vc[TRC + t][k][j][i] = 1.0;
                        }
                        else
                        {
                            d->Vc[TRC + t][k][j][i] = 0.0;
                        }
                    }
                }
            }
            else
            {  /* -- reflective boundary for r > theta_rad --*/
                d->Vc[RHO][k][j][i] =  d->Vc[RHO][k][j][2*IBEG - i - 1];
                d->Vc[VX1][k][j][i] =  d->Vc[VX1][k][j][2*IBEG - i - 1];
                d->Vc[VX2][k][j][i] = -d->Vc[VX2][k][j][2*IBEG - i - 1];
                d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][j][2*IBEG - i - 1];
                if (NTRACER > 0)
                {
                    d->Vc[TRC][k][j][i] =  d->Vc[TRC][k][j][2*IBEG - i - 1]; /* tracer particle over jet */
                }
            }

        }

    }
}

#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 *********************************************************************** */
{
    
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the graviational potential as function of the coordinates.
 *
 *********************************************************************** */

{
   // double Phi_0;
   // Phi_0 = -1.5*g_inputParam[B_EXPONENT]/g_gamma;
   // return Phi_0*log(1+(x1/g_inputParam[R_CORE])*(x1/g_inputParam[R_CORE]));
    double rho, g_gamma, r;
    g_gamma = 5.0/3.0;
    r = x1;
    rho = pow(1 + pow(r/g_inputParam[R_CORE],2),-1.5*g_inputParam[B_EXPONENT]);
    return -log(rho)/g_gamma;
    
}
#endif
