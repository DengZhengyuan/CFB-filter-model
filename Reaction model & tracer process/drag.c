// -------------------------------------------------------------
//                      Gidaspow drag + Hd
// -------------------------------------------------------------
DEFINE_EXCHANGE_PROPERTY(Gidaspow_Hd_drag, cell, thread_mix, s_phase, f_phase)
{
    Thread *thread_g, *thread_s;

    real rho_g, rho_s, mu_g, epsilon_g, epsilon_s, cell_vol, u_slip;
    real C_D, num_Re, k_g_s_HG, coef_Hd, k_g_s;

    thread_g = THREAD_SUB_THREAD(thread_mix, s_phase); // gas phase
    thread_s = THREAD_SUB_THREAD(thread_mix, f_phase); // solid phase

    rho_g = C_R(cell, thread_g);           // density of gas
    rho_s = C_R(cell, thread_s);           // density of solid
    mu_g = C_MU_L(cell, thread_g);         // viscosity of gas
    epsilon_g = C_VOF(cell, thread_g);     // voidage
    epsilon_s = 1. - epsilon_g + 1.e-8;    // solids holdup
    cell_vol = C_VOLUME(cell, thread_mix); // cell volume
    u_slip = C_UDMI(cell, thread_g, 0);    // slip velocity

    num_Re = rho_g * u_slip * diam_p / mu_g; // Reynolds number
    // drag coefficient
    C_D = 24. * (1. + 0.15 * pow(epsilon_g * num_Re, 0.687)) / (epsilon_g * num_Re + 1.e-8); 

    if (epsilon_g > 0.8)
    { // Wen-Yu
        k_g_s_HG = 0.75 * epsilon_g * epsilon_s * rho_g * u_slip * C_D * pow(epsilon_g, -2.65) / diam_p;
    }
    else
    { // Ergun
        k_g_s_HG = 150. * (mu_g * pow(epsilon_s, 2.) / (epsilon_g * pow(diam_p, 2.))) + 1.75 * u_slip * (rho_g * epsilon_s) / diam_p;
    }

    coef_Hd = 1.;
    k_g_s = coef_Hd * k_g_s_HG;

    C_UDMI(cell, thread_mix, 4) = coef_Hd;

    return k_g_s;
}