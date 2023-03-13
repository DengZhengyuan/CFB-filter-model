#define oz_conc 1.6653e-4  /* initial ozone concentration, 100 ppm -> %mass */
#define adjust_ID 4

#include "udf.h"

DEFINE_ON_DEMAND(test_ML) 
{
	Domain *domain, *domain_gas, *domain_solids; 
    Thread *thread;
    face_t face;
    cell_t cell;
    
    domain = Get_Domain(1); // mixture domain
	domain_gas = Get_Domain(2); // gas domain
	domain_solids = Get_Domain(3); // solids domain

    /* Loop over all cell threads in domain */
	thread_loop_c(thread, domain)
	{ /* Loop over all cells */
		begin_c_loop(cell, thread)
		{ 
			// solids holdup gradient []
			C_UDMI(cell, thread, 5) = C_UDSI_G(cell, thread, 1)[0];
			// C_UDMI(cell, thread, 6) = C_UDSI_G(cell, thread, 1)[1];
			// test slip velocity gradient []
			// C_UDMI(cell, thread, 7) = C_UDSI_G(cell, thread, 2)[0];
			// C_UDMI(cell, thread, 8) = C_UDSI_G(cell, thread, 2)[1];
			// pressure gradient []
			// C_UDMI(cell, thread, 9) = C_P_G(cell, thread)[0];
			// C_UDMI(cell, thread, 10) = C_P_G(cell, thread)[1];
			// test ozone gradient [x]
			// C_UDMI(cell, thread, 11) = C_UDSI_G(cell, thread, 0)[0] / oz_conc;
			// C_UDMI(cell, thread, 12) = C_UDSI_G(cell, thread, 0)[1] / oz_conc;
		}
		end_c_loop(cell, thread);
	}
	printf("On demand executed. \n");
}

DEFINE_ON_DEMAND(thread_cell_loop)
{
	Domain *domain;
	Thread *thread, *thread_gas, *thread_solids;
	face_t face;
	cell_t cell;

	domain = Get_Domain(1); // mixture domain
	thread = Lookup_Thread(domain, adjust_ID);
	thread_solids = THREAD_SUB_THREAD(thread, 1);
	thread_gas = THREAD_SUB_THREAD(thread, 0);

	/* Loop over all cells */
	begin_c_loop(cell, thread)
	{
		// solids holdup gradient []
		// C_UDMI(cell, thread, 5) = C_UDSI_G(cell, thread, 1)[0];
		C_UDMI(cell, thread, 6) = C_UDSI_G(cell, thread, 1)[0];
		// test slip velocity gradient []
		C_UDMI(cell, thread, 7) = C_UDSI_G(cell, thread, 2)[0];
		C_UDMI(cell, thread, 8) = C_UDSI_G(cell, thread, 2)[1];
		// pressure gradient []
		C_UDMI(cell, thread, 9) = C_P_G(cell, thread)[0];
		C_UDMI(cell, thread, 10) = C_P_G(cell, thread)[1];
		// test ozone gradient [x]
		C_UDMI(cell, thread, 11) = (C_UDSI_G(cell, thread_gas, 0)[0]) / oz_conc;
		C_UDMI(cell, thread, 12) = C_UDSI_G(cell, thread, 0)[1] / oz_conc;
	}
	end_c_loop(cell, thread)
	printf("On demand executed. \n");
}
