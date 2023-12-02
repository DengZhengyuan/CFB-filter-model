#define gt_conc 1.6653e-4 // initial gtone concentration, 100 ppm -> %mass
#define domain_ID 4
#define uds_tracer_ID_gas 3
#define uds_tracer_ID_solids 4

#include "udf.h"

DEFINE_EXECUTE_AT_END(write_tracer)
{
    if (N_TIME % 10 == 0)
    {
#if !RP_HOST
        Domain *domain;
        Thread *thread, *thread_gas, *thread_solids;
        face_t face;
        cell_t cell;
#else
        int i;
#endif

#if !RP_NODE
        FILE *fp = NULL;
        char filename[] = "tracer_out.csv"; // output file name
#endif

        int size;          // data passing variables
        real *array_gt;    // tracer in gas phase
        real *array_st;    // tracer in solids phase
        real *array_pos_x; // x-coordinate
        real *array_pos_y; // y-coordinate
        real *array_es;    // solids holdup
        real *array_Ug;    // gas phase velocity in y-direction
        real *array_Us;    // solids phase velocity in y-direction
        real x[ND_ND];     // get coordinates
        int pe;

#if !RP_HOST
        domain = Get_Domain(1); // mixture domain
        thread = Lookup_Thread(domain, domain_ID);
        thread_solids = THREAD_SUB_THREAD(thread, 1);
        thread_gas = THREAD_SUB_THREAD(thread, 0);
#endif

#if !RP_NODE
        fp = fopen(filename, "a");
        Message("\nWriting data to %s...", filename);
#endif

#if RP_NODE
        // Each Node loads up its data passing array_gt
        size = THREAD_N_ELEMENTS_INT(thread);
        array_gt = (real *)malloc(size * sizeof(real));
        array_st = (real *)malloc(size * sizeof(real));
        array_pos_x = (real *)malloc(size * sizeof(real));
        array_pos_y = (real *)malloc(size * sizeof(real));
        array_es = (real *)malloc(size * sizeof(real));
        array_Ug = (real *)malloc(size * sizeof(real));
        array_Us = (real *)malloc(size * sizeof(real));
        begin_c_loop_int(cell, thread_gas)
        {
            array_gt[cell] = C_UDSI(cell, thread_gas, uds_tracer_ID_gas); // mass fraction of tracer
            array_es[cell] = 1. - C_VOF(cell, thread_gas);                // solids holdup
            C_CENTROID(x, cell, thread_gas);                              // get coordinates
            array_pos_x[cell] = x[0];                                     // x-coordinate
            array_pos_y[cell] = x[1];                                     // y-coordinate
            array_Ug[cell] = C_V(cell, thread_gas);                       // gas velocity in y-direction
        }
        end_c_loop_int(cell, thread_gas)

            begin_c_loop_int(cell, thread_solids)
        {
            array_st[cell] = C_UDSI(cell, thread_solids, uds_tracer_ID_solids); // mass fraction of tracer
            array_Us[cell] = C_V(cell, thread_solids);                          // solids velocity in y-direction
        }
        end_c_loop_int(cell, thread_solids)

            // Set pe to destination node
            // If on node_0 send data to host
            // Else send to node_0 because
            // compute nodes connect to node_0 & node_0 to host
            pe = (I_AM_NODE_ZERO_P) ? node_host : node_zero;
        PRF_CSEND_INT(pe, &size, 1, myid);
        PRF_CSEND_REAL(pe, array_gt, size, myid);
        PRF_CSEND_REAL(pe, array_st, size, myid);
        PRF_CSEND_REAL(pe, array_es, size, myid);
        PRF_CSEND_REAL(pe, array_pos_x, size, myid);
        PRF_CSEND_REAL(pe, array_pos_y, size, myid);
        PRF_CSEND_REAL(pe, array_Ug, size, myid);
        PRF_CSEND_REAL(pe, array_Us, size, myid);
        free(array_gt);
        free(array_st);
        free(array_es);
        free(array_pos_x);
        free(array_pos_y);
        free(array_Ug);
        free(array_Us);
        // free array_gt on nodes after data sent
        // node_0 now collect data sent by other compute nodes
        // and sends it straight on to the host
        if (I_AM_NODE_ZERO_P)
            compute_node_loop_not_zero(pe)
            {
                PRF_CRECV_INT(pe, &size, 1, pe);
                array_gt = (real *)malloc(size * sizeof(real));
                array_st = (real *)malloc(size * sizeof(real));
                array_es = (real *)malloc(size * sizeof(real));
                array_pos_x = (real *)malloc(size * sizeof(real));
                array_pos_y = (real *)malloc(size * sizeof(real));
                array_Ug = (real *)malloc(size * sizeof(real));
                array_Us = (real *)malloc(size * sizeof(real));
                PRF_CRECV_REAL(pe, array_gt, size, pe);
                PRF_CRECV_REAL(pe, array_st, size, pe);
                PRF_CRECV_REAL(pe, array_es, size, pe);
                PRF_CRECV_REAL(pe, array_pos_x, size, pe);
                PRF_CRECV_REAL(pe, array_pos_y, size, pe);
                PRF_CRECV_REAL(pe, array_Ug, size, pe);
                PRF_CRECV_REAL(pe, array_Us, size, pe);
                PRF_CSEND_INT(node_host, &size, 1, myid);
                PRF_CSEND_REAL(node_host, array_gt, size, myid);
                PRF_CSEND_REAL(node_host, array_st, size, myid);
                PRF_CSEND_REAL(node_host, array_es, size, myid);
                PRF_CSEND_REAL(node_host, array_pos_x, size, myid);
                PRF_CSEND_REAL(node_host, array_pos_y, size, myid);
                PRF_CSEND_REAL(node_host, array_Ug, size, myid);
                PRF_CSEND_REAL(node_host, array_Us, size, myid);
                free((char *)array_gt);
                free((char *)array_st);
                free((char *)array_es);
                free((char *)array_pos_x);
                free((char *)array_pos_y);
                free((char *)array_Ug);
                free((char *)array_Us);
            }
#endif // RP_NODE

#if RP_HOST
        compute_node_loop(pe) // only acts as a counter in this loop
        {
            // Receive data sent by each node and write it out to the file
            PRF_CRECV_INT(node_zero, &size, 1, node_zero);
            array_gt = (real *)malloc(size * sizeof(real));
            array_st = (real *)malloc(size * sizeof(real));
            array_es = (real *)malloc(size * sizeof(real));
            array_pos_x = (real *)malloc(size * sizeof(real));
            array_pos_y = (real *)malloc(size * sizeof(real));
            array_Ug = (real *)malloc(size * sizeof(real));
            array_Us = (real *)malloc(size * sizeof(real));
            PRF_CRECV_REAL(node_zero, array_gt, size, node_zero);
            PRF_CRECV_REAL(node_zero, array_st, size, node_zero);
            PRF_CRECV_REAL(node_zero, array_es, size, node_zero);
            PRF_CRECV_REAL(node_zero, array_pos_x, size, node_zero);
            PRF_CRECV_REAL(node_zero, array_pos_y, size, node_zero);
            PRF_CRECV_REAL(node_zero, array_Ug, size, node_zero);
            PRF_CRECV_REAL(node_zero, array_Us, size, node_zero);
            for (i = 0; i < size; i++)
            {
                if (array_pos_y[i] > 9.999 && array_pos_y[i] < 10.001) // 10 meters, outlet
                    fprintf(fp, "%g, %g, %g, %g, %g, %g, %g, %g\n", CURRENT_TIME, array_pos_x[i], array_pos_y[i], array_es[i], array_Ug[i], array_Us[i], array_gt[i], array_st[i]);
                if (array_pos_y[i] > 4.999 && array_pos_y[i] < 5.001) // 5 meters
                    fprintf(fp, "%g, %g, %g, %g, %g, %g, %g, %g\n", CURRENT_TIME, array_pos_x[i], array_pos_y[i], array_es[i], array_Ug[i], array_Us[i], array_gt[i], array_st[i]);
                if (array_pos_y[i] > 1.999 && array_pos_y[i] < 2.001) // 2 meters
                    fprintf(fp, "%g, %g, %g, %g, %g, %g, %g, %g\n", CURRENT_TIME, array_pos_x[i], array_pos_y[i], array_es[i], array_Ug[i], array_Us[i], array_gt[i], array_st[i]);
            }
            free(array_gt);
            free(array_st);
            free(array_es);
            free(array_pos_x);
            free(array_pos_y);
            free(array_Ug);
            free(array_Us);
        }
#endif // RP_HOST

#if !RP_NODE
        fclose(fp);
        Message("Done\n");
#endif
    }
}
