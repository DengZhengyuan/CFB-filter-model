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

#define diam_p 70.e-6      /* particle diameter [m] */
#define oz_conc 1.6653e-4  /* initial ozone concentration, 100 ppm -> %mass */
#define k_v 49.2 /* reaction constant based on catalysts volume, [/s] */
#define D_ozone 1.48535e-5 /* mass diffusivity of ozone in air [m^2/s] */
#define density_gas 1.225        /* density of air, [kg/m^3] */
#define density_solids 1780.     /* density of solids, [kg/m^3] */
#define adjust_ID 491

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
#include <sys/socket.h>
#include <sys/un.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <zconf.h>
#include <vector.c>

const u_int8_t PACKETSEP = '\0';
const u_int8_t CMD_KILL = 0;	// kill server
const u_int8_t CMD_FEATURE = 1; // Calculate Feature Data
const u_int8_t CMD_RESULT = 2;  // Return Value

#define DBG_ERROR(...)                                          \
  fprintf(stderr, "Error Log At: %s %d\n", __FILE__, __LINE__); \
  fprintf(stderr, __VA_ARGS__);                                 \
  fprintf(stderr, "\n");
#define DBG_LOG(...)

/**
 * Open a socket connection
 * @param socket_path the path of us
 * @return File discriptor, -1 if failed
 */
int openSocket(char *socket_path)
{
	struct sockaddr_un addr;
	int fd;
	/* Create Domain Socket */
	if ((fd = socket(AF_UNIX, SOCK_STREAM, 0)) == -1)
	{
		DBG_ERROR("socket connection error %d", fd);
		return -1;
	}

	/* setup socket */
	DBG_LOG("setup socket");
	memset(&addr, 0, sizeof(addr));
	DBG_LOG("memset");
	addr.sun_family = AF_UNIX;
	strncpy(addr.sun_path, socket_path, sizeof(addr.sun_path) - 1);
	DBG_LOG("strncpy");
	if (connect(fd, (struct sockaddr *)&addr, sizeof(addr)) == -1)
	{
		DBG_ERROR("socket connection error %d", -1);
		return -1;
	}
	DBG_LOG("Here");
	return fd;
}

/**
 * Read a whole packet from socket data
 * This function will read the whole packet and return an array
 * allocated using malloc. Don't forget to call free when finished.
 *
 * @param fd
 * @param bif
 * @param bufLen
 * @param result dataLen or -1 when error occurred
 */
int readPacket(int fd, u_int8_t *command, u_int32_t *rltDataLen, void **rltData)
{
  /*1.Read till a PACKETSEP*/
  while (1)
  {
    u_int8_t firstByte;
    DBG_LOG("Read first Byte");
    ssize_t rc = recv(fd, &firstByte, 1, MSG_WAITALL); /* Read first byte */
    if (rc != 1)
    {
      DBG_ERROR("HEADER Parse ERROR %d", rc);
      return -1;
    }
    else if (firstByte == PACKETSEP)
      break;
  }
  /* 2. Read command */
  ssize_t rc = recv(fd, command, 1, MSG_WAITALL);
  DBG_LOG("Read command %d", *command);

  /* 3. Read packet length */
  rc = recv(fd, rltDataLen, 4, MSG_WAITALL);
  DBG_LOG("Read packet length %d", *rltDataLen);

  if (rc != 4)
  {
    DBG_ERROR("RECV ERROR %d", rc);
    return -1;
  }

  /* 4. Read Data */
  *rltData = (char *)malloc(*rltDataLen);
  rc = recv(fd, *rltData, *rltDataLen, MSG_WAITALL);
  DBG_LOG("Read data %d", *rltDataLen);

  if (rc != *rltDataLen)
  {
    DBG_ERROR("RECV ERROR %d", rc);
    return -1;
  }
  return *rltDataLen;
}

/**
 * Write data to socket
 * Packet format is: PACKETSEP(8b) command(8b) packetLen(32b) packetData(packetLen b)
 * Packet header is automatically added
 * @param fd Socket file discriptor
 * @param buf Data to be sent
 * @param bufLen Length of buf
 * @return length of sent bytes, -1 if failed
 */
ssize_t writeSocket(int fd, u_int8_t command, u_int32_t bufLen, void *buf)
{
  write(fd, &PACKETSEP, 1); // write header
  write(fd, &command, 1);   // write command
  write(fd, &bufLen, 4);    // write bufferLen

  ssize_t rc = write(fd, buf, bufLen); // write buffer

  if (rc == -1)
  {
    DBG_ERROR("write error return code is -1");
    return -1;
  }
  else if (rc != bufLen)
  {
    DBG_ERROR("partial write %d", rc);
    return -1;
  }
  /* write success */
  return bufLen;
}

int closeSocket(int fd)
{
	return close(fd);
}

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
#include "udf.h"
#include "stdio.h"
#include "math.h"
#include "sg_pb.h"
#include "sg_mphase.h"
#include "mem.h"
#include "vector.h"


int fd = -1;

int isFirstRun = 1;
Vector vec_nt,                             // 1
    vec_Es, vec_Es_dx, vec_Es_dy,          // 2
    vec_Uslip, vec_Uslip_dx, vec_Uslip_dy, // 3
    vec_P_dx, vec_P_dy,                    // 4
    vec_Oz, vec_Oz_dx, vec_Oz_dy,          // 6
    vec_beta;                              // 7

double *allData = NULL;
double *hdData = NULL;
double *hrData = NULL;

int i;
u_int32_t allDataLen = 0;
u_int32_t hdDataLen = 0;
u_int32_t hrDataLen = 0;


DEFINE_EXECUTE_AT_END(feature_loop)
{
  // printf("Start feature loop.\n ");

  if (myid == node_host)
  {
    DBG_LOG("I'm host and I shouldn't do calculation!");
    return;
  }

  if (isFirstRun)
  {
    DBG_LOG("First Run, H_d = 1. \n");
    // printf("First Run. \n");
    isFirstRun = 0;
    vector_setup(&vec_nt, 10, sizeof(real)); // 1
    vector_setup(&vec_Es, 10, sizeof(real));
    vector_setup(&vec_Es_dx, 10, sizeof(real)); 
    vector_setup(&vec_Es_dy, 10, sizeof(real)); // 2
    vector_setup(&vec_Uslip, 10, sizeof(real));
    vector_setup(&vec_Uslip_dx, 10, sizeof(real));
    vector_setup(&vec_Uslip_dy, 10, sizeof(real)); // 3
    vector_setup(&vec_P_dx, 10, sizeof(real));
    vector_setup(&vec_P_dy, 10, sizeof(real)); //4 
    vector_setup(&vec_Oz, 10, sizeof(real));
    vector_setup(&vec_Oz_dx, 10, sizeof(real));
    vector_setup(&vec_Oz_dy, 10, sizeof(real)); // 6
    vector_setup(&vec_beta, 10, sizeof(real)); // 7
  }
  else
  {
    free(allData);
    free(hdData);
    free(hrData);
    allData = NULL;
    allDataLen = -1;
    hdData = NULL;
    hdDataLen = -1;
    hrData = NULL;
    hrDataLen = -1;
    //
    vector_clear(&vec_nt);
    vector_clear(&vec_Es);
    vector_clear(&vec_Es_dx);
    vector_clear(&vec_Es_dy);
    vector_clear(&vec_Uslip);
    vector_clear(&vec_Uslip_dx);
    vector_clear(&vec_Uslip_dy);
    vector_clear(&vec_P_dx);
    vector_clear(&vec_P_dy);
    vector_clear(&vec_Oz);
    vector_clear(&vec_Oz_dx);
    vector_clear(&vec_Oz_dy);
    vector_clear(&vec_beta);
  }

  /* Define variables */
  real filter_size,
      es, es_dx, es_dy,
      Uslip, Uslip_dx, Uslip_dy,
      p_dx, p_dy,
      oz, oz_dx, oz_dy, 
      beta;
  // real sum_nt = 0.;

  Domain *domain, *domain_gas, *domain_solids; 
	Thread *thread, *thread_gas, *thread_solids; 
  // face_t face;
  cell_t cell;

	domain = Get_Domain(1); // mixture domain
	domain_gas = Get_Domain(2); // gas domain
	domain_solids = Get_Domain(3); // solids domain
 
	if(!domain_gas)
		DBG_LOG("Domain gas failed");
	if(!domain_solids)
		DBG_LOG("Domain solids failed");

  thread = Lookup_Thread(domain, adjust_ID);
	thread_gas = Lookup_Thread(domain_gas, adjust_ID); /*  gas phase  */
	thread_solids = Lookup_Thread(domain_solids, adjust_ID); /*  solid phase */

  if (!thread_gas)
    DBG_LOG("Gas thread failed");
  if (!thread_solids)
    DBG_LOG("Solids thread failed");

  if (fd == -1)
	{
		// printf("OpenSocket.\n");
		DBG_LOG("Before open Socket");
		fd = openSocket((char *)"/tmp/welfkewgsocket");
	}
	DBG_LOG("After open socket");

  DBG_LOG("Begin Loop");
  // printf("Begin Loop. \n");

  begin_c_loop(cell, thread)
  {
    // filter size (nt)
    filter_size = pow(C_VOLUME(cell, thread_gas), 1./3.) / diam_p;
    // sum_nt += filter_size;
    vector_push_back(&vec_nt, &filter_size);

    // slip velocity (Uslip)
    Uslip = C_UDMI(cell, thread, 0);
    vector_push_back(&vec_Uslip, &Uslip);
    // gradient
    Uslip_dx = C_UDSI_G(cell, thread, 2)[0] * C_UDSI_G(cell, thread, 2)[2];
    vector_push_back(&vec_Uslip_dx, &Uslip_dx);
    Uslip_dy = C_UDSI_G(cell, thread, 2)[1];
    vector_push_back(&vec_Uslip_dy, &Uslip_dy);

    // pressure gradient (p_dx, p_dy)
    p_dx = C_P_G(cell, thread)[0] * C_P_G(cell, thread)[2];
    vector_push_back(&vec_P_dx, &p_dx);
    p_dy = C_P_G(cell, thread)[1];
    vector_push_back(&vec_P_dy, &p_dy);

    // ozone concentration (oz)
    oz = C_UDSI(cell, thread_gas, 0) / oz_conc;
    vector_push_back(&vec_Oz, &oz);
    // ozone concentration gradient (Oz_dx, Oz_dy)
    oz_dx = (C_UDSI_G(cell, thread_gas, 0)[0]*C_UDSI_G(cell, thread_gas, 0)[2]) / oz_conc;
    vector_push_back(&vec_Oz_dx, &oz_dx);
    oz_dy = C_UDSI_G(cell, thread_gas, 0)[1] / oz_conc;
    vector_push_back(&vec_Oz_dy, &oz_dy);

    // solids holdup (Es)
    es = 1.- C_VOF(cell, thread_gas);
    vector_push_back(&vec_Es, &es);
    // solids holdup gradient (Es_dx, Es_dy)
    es_dx = C_UDSI_G(cell, thread, 1)[0] * C_UDSI_G(cell, thread, 1)[2];
    vector_push_back(&vec_Es_dx, &es_dx);
    es_dy = C_UDSI_G(cell, thread, 1)[1];
    vector_push_back(&vec_Es_dy, &es_dy);

    // drag coefficient 
    beta = C_UDMI(cell, thread, 5);
    vector_push_back(&vec_beta, &beta);
  }
  end_c_loop(cell, thread)

	// printf("End Loop. \n");
	DBG_LOG("End Loop");
	DBG_LOG("");

	/*sending data to server */
	int32_t sizeCache = vec_nt.size;
  // printf("The vector size is %d. \n", vec_nt.size);
	DBG_LOG("Writing arrLen====================%d", sizeCache);
	writeSocket(fd, CMD_FEATURE, 4, &sizeCache);
  DBG_LOG("Writing parameters %d", 8 * vec_nt.size);
	writeSocket(fd, CMD_FEATURE, 8 * vec_nt.size, vec_nt.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_Es.size, vec_Es.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_Es_dx.size, vec_Es_dx.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_Es_dy.size, vec_Es_dy.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_Uslip.size, vec_Uslip.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_Uslip_dx.size, vec_Uslip_dx.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_Uslip_dy.size, vec_Uslip_dy.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_P_dx.size, vec_P_dx.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_P_dy.size, vec_P_dy.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_Oz.size, vec_Oz.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_Oz_dx.size, vec_Oz_dx.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_Oz_dy.size, vec_Oz_dy.data);
	writeSocket(fd, CMD_FEATURE, 8 * vec_beta.size, vec_beta.data);

	/*reading results*/
  DBG_LOG("Read Packet");
	u_int8_t command = 0;
	readPacket(fd, &command, &hdDataLen, &hdData);
	readPacket(fd, &command, &hrDataLen, &hrData);
	DBG_LOG("Reading complete");

}

DEFINE_EXECUTE_AT_EXIT(destroy_vector){
  vector_destroy(&vec_nt);
  vector_destroy(&vec_Es);
  vector_destroy(&vec_Es_dx);
  vector_destroy(&vec_Es_dy);
  vector_destroy(&vec_Uslip);
  vector_destroy(&vec_Uslip_dx);
  vector_destroy(&vec_Uslip_dy);
  vector_destroy(&vec_P_dx);
  vector_destroy(&vec_P_dy);
  vector_destroy(&vec_Oz);
  vector_destroy(&vec_Oz_dx);
  vector_destroy(&vec_Oz_dy);
  vector_destroy(&vec_beta);
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

    real coord_x[ND_ND];
    C_CENTROID(coord_x, cell, thread_mix);
    
    thread_g = THREAD_SUB_THREAD(thread_mix, s_phase); /* gas phase */
    thread_s = THREAD_SUB_THREAD(thread_mix, f_phase); /* solid phase */

    rho_g = C_R(cell, thread_g);             /* density of gas */
    rho_s = C_R(cell, thread_s);             /* density of solid */
    mu_g = C_MU_L(cell, thread_g);           /* viscosity of gas */
    epsilon_g = C_VOF(cell, thread_g);       /* voidage */
    epsilon_s = 1. - epsilon_g + 1.e-8; /* solids holdup */
    cell_vol = C_VOLUME(cell, thread_mix);   /* cell volume */
    u_slip = C_UDMI(cell, thread_g, 0);      /* slip velocity */

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
    if (hdData)
    {
      coef_Hd = hdData[cell];

      if (coef_Hd < 0.03)
        coef_Hd = 0.03;

      if (coef_Hd > 1.)
        coef_Hd = 1.;

      // if (epsilon_s >= 0.54)
      //   coef_Hd = 1.;
    }
    else
    {
      coef_Hd = 1.;
    }

    // if (coord_x[0] > 0.0381 || coord_x[1] < 0.05)
    // {
    //   coef_Hd = 0.8;
    // }

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
    
    /* --------------------------------- */
    if (hrData)
    {
      coef_Hr = hrData[cell];
      if (coef_Hr < 0.03)
        coef_Hr = 0.03;

      if (coef_Hr > 1.)
        coef_Hr = 1.;
    }
    else
    {
      coef_Hr = 1.;
    }
    /*---- calculate the reaction rate ----*/
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
