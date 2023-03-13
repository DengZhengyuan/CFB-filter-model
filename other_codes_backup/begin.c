/*---------------------------------------------------------------------------
    Created by Zhengyuan (Zach) Deng
    Modified on July 15, 2022
    
    This UDF contains 
        - 
---------------------------------------------------------------------------*/

/*
-------------------------------------------------------

    UDS_0 -> %mass of ozone in gas phase
    UDS_1 -> solids holdup in mixture phase
    UDS_2 -> slip velocity in mixture phase

-------------------------------------------------------

    UDM_0 -> slip velocity, m/s
    UDM_1 -> reaction rate, kg/m^3/s
    UDM_2 -> ozone concentration, kg/m^3/s
--------------
    UDM_3 -> coefficient of Hr
    UDM_4 -> coefficient of Hd

-------------------------------------------------------
*/
// #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
// #pragma GCC diagnostic ignored "-Wformat="
// #pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
// #pragma GCC diagnostic ignored "-Wunused-variable"

#define ID_outlet 24             /* zone ID of riser's outlet */
#define e_s_0 0.5                /* solids holdup at inlet */
#define A_riser 0.004545       /* cross-sectional area of riser, [m^2] */
#define A_inlet_solids 0.00201655 /* area of solids inlet, [m^2] */

#define diam_p 70.e-6      /* particle diameter [m] */
#define oz_conc 1.6653e-4  /* initial ozone concentration, 100 ppm -> %mass */
#define k_v 49.2 /* reaction constant based on catalysts volume, [/s] */
#define D_ozone 1.48535e-5 /* mass diffusivity of ozone in air [m^2/s] */
#define density_gas 1.225        /* density of air, [kg/m^3] */
#define density_solids 1780.     /* density of solids, [kg/m^3] */
#define diam_riser 0.0762 // diameter of riser, [m]

#include "udf.h"
#include "stdio.h"
#include "math.h"

/*---------------------------------------------------------------------------*/
/*-------------- set velocity inlet of solids at solids inlet ---------------*/
/*---------------------------------------------------------------------------*/
DEFINE_PROFILE(Us_at_solids_inlet, thread, i)
{
#if !RP_HOST
    Domain *domain;
    face_t face;
    cell_t cell;
    Thread *thread_mix, *thread_gas, *thread_solids;

    real flux = 0.0;
    real vfr_s_outlet = 0.0;

    domain = Get_Domain(1);
    thread_mix = Lookup_Thread(domain, ID_outlet);
    thread_gas = THREAD_SUB_THREAD(thread_mix, 0);
    thread_solids = THREAD_SUB_THREAD(thread_mix, 1);
#endif

#if !RP_HOST
    begin_f_loop(face, thread_solids)
    {
        flux += F_FLUX(face, thread_solids);
    }
    end_f_loop(face, thread_solids)

        flux = PRF_GRSUM1(flux);
    vfr_s_outlet = flux / density_solids;

    begin_f_loop(face, thread)
    {
        F_PROFILE(face, thread, i) = vfr_s_outlet / (A_inlet_solids * e_s_0);
    }
    end_f_loop(face, thread)
#endif
}

/* ------------------------------------------------------------- */
/*                      Gidaspow drag + Hd                       */
/* ------------------------------------------------------------- */
DEFINE_EXCHANGE_PROPERTY(Gidaspow_Hd_drag, cell, thread_mix, s_phase, f_phase)
{
    Thread *thread_g, *thread_s;

    real rho_g,
        rho_s,
        mu_g,
        epsilon_g,
        epsilon_s,
        cell_vol,
        u_slip;
    real C_D,
        num_Re,
        k_g_s_HG,
        coef_Hd,
        k_g_s;

    thread_g = THREAD_SUB_THREAD(thread_mix, s_phase); /* gas phase */
    thread_s = THREAD_SUB_THREAD(thread_mix, f_phase); /* solid phase */

    rho_g = C_R(cell, thread_g);           /* density of gas */
    rho_s = C_R(cell, thread_s);           /* density of solid */
    mu_g = C_MU_L(cell, thread_g);         /* viscosity of gas */
    epsilon_g = C_VOF(cell, thread_g);     /* voidage */
    epsilon_s = 1. - epsilon_g + 1.e-8;    /* solids holdup */
    cell_vol = C_VOLUME(cell, thread_mix); /* cell volume */
    u_slip = C_UDMI(cell, thread_g, 0);    /* slip velocity */

    num_Re = rho_g * u_slip * diam_p / mu_g; /* Reynolds number */
    C_D = 24. * (1. + 0.15 * pow(epsilon_g * num_Re, 0.687)) / (epsilon_g * num_Re + 1.e-8); /* drag coefficient */

    if (epsilon_g > 0.8)
    { /* Wen-Yu */
      k_g_s_HG = 0.75 * epsilon_g * epsilon_s * rho_g * u_slip * C_D * pow(epsilon_g, -2.65) / diam_p;
    }
    else
    { /* Ergun */
      k_g_s_HG = 150. * (mu_g * pow(epsilon_s, 2.) / (epsilon_g * pow(diam_p, 2.))) + 1.75 * u_slip * (rho_g * epsilon_s) / diam_p;
    }

    /* --------------------------------- */
    // if (hdData)
    // {
    //   coef_Hd = hdData[cell];

    //   if (coef_Hd < 0.03)
    //     coef_Hd = 0.03;

    //   if (coef_Hd > 1.5)
    //     coef_Hd = 1.5;

    //   if (epsilon_g >= 0.995)
    //     coef_Hd = 1.;
    // }
    // else
    // {
    //   coef_Hd = 1.;
    // }
    coef_Hd = 1.;
    k_g_s = coef_Hd * k_g_s_HG;

    C_UDMI(cell, thread_mix, 4) = coef_Hd;
    C_UDMI(cell, thread_mix, 5) = k_g_s_HG;

    return k_g_s;
}

/* ------------------------------------------------------------- */
/*                    reaction term of ozone                     */
/* ------------------------------------------------------------- */
DEFINE_SOURCE(rxn_ozone, cell, thread, dS, eqn)
{
    real rxn_rate;
    real Y_ozone;
    real rho_g;
    real epsilon_g, epsilon_s;
    real coef_Hr;
    real k;

    Y_ozone = C_UDSI(cell, thread, 0); /* mole of ozone in the cell */
    epsilon_g = C_VOF(cell, thread);   /* volume fraction of gas */
    epsilon_s = 1. - epsilon_g;        /* volume fraction of solids */
    rho_g = C_R(cell, thread);         /* gas density */
    
    /*---- calculate the reaction rate ----*/
    coef_Hr = 1.;
    k = coef_Hr * k_v;
    rxn_rate = -k * rho_g * epsilon_s * Y_ozone;
    dS[eqn] = -k * rho_g * epsilon_s;
    C_UDMI(cell, thread, 1) = rxn_rate;
    C_UDMI(cell, thread, 2) = Y_ozone/oz_conc;
    C_UDMI(cell, thread, 3) = coef_Hr;

    return rxn_rate;
}

/* ------------------------------------------------------------- */
/*                     diffusivity of ozone                      */
/* ------------------------------------------------------------- */
DEFINE_DIFFUSIVITY(diff_ozone_laminar, cell, thread, i)
{
    Thread *thread_mix, *thread_solids, *thread_gas;

    real diff, D_AB;
    real rho_g, u_slip;
    real epsilon_g, epsilon_s;
    real Y_oz_g;

    // real coord_x[ND_ND];
    // C_CENTROID(coord_x, cell, thread_mix);

    /*---------------------------------*/
    /*--------- slip velocity ---------*/
    /* find the solids thread */
    thread_mix = THREAD_SUPER_THREAD(thread);
    thread_solids = THREAD_SUB_THREAD(thread_mix, 1);
    thread_gas = THREAD_SUB_THREAD(thread_mix, 0);

#if RP_2D
    u_slip = pow(
        pow(C_U(cell, thread) - C_U(cell, thread_solids), 2.) +
            pow(C_V(cell, thread) - C_V(cell, thread_solids), 2.),
        .5);
#endif

#if RP_3D
    u_slip = pow(
        pow(C_U(cell, thread) - C_U(cell, thread_solids), 2.) +
            pow(C_V(cell, thread) - C_V(cell, thread_solids), 2.) +
            pow(C_W(cell, thread) - C_W(cell, thread_solids), 2.),
        .5);
#endif
    /* UDM_0 -> slip velocity */
    C_UDMI(cell, thread, 0) = u_slip;

    /*--------- slip velocity ---------*/
    /*---------------------------------*/

    /* volume fraction of solids and gas */
    epsilon_g = C_VOF(cell, thread);
    epsilon_s = 1. - epsilon_g;
    Y_oz_g = C_UDSI(cell, thread_gas, 0); /* %mass of ozone in the gas phase */

    /*---------------------------------*/
    // if (sqrt((x[0])^2.+(x[2])^2.) > diam_riser/2)
    // {
    //     C_UDSI(cell, thread_mix, 1) = 0.5;
    //     C_UDSI(cell, thread_mix, 2) = 0.5;
    // }
    // else{
    //     C_UDSI(cell, thread_mix, 1) = epsilon_s;
    //     C_UDSI(cell, thread_mix, 2) = u_slip;
    // }
    C_UDSI(cell, thread_mix, 1) = epsilon_s;
    C_UDSI(cell, thread_mix, 2) = u_slip;
    /*---------------------------------*/

    /* define some parameters and properties */
    rho_g = C_R(cell, thread);
    D_AB = D_ozone * (1. - pow(epsilon_s, 0.5)) / epsilon_g;

    if (Y_oz_g < (0.01 * oz_conc))
    {
        diff = 0.;
    }
    else
    {
        diff = rho_g * D_AB;
    }

    return diff;
}

