/* Start C_IO_TIMER                                                 */
/*                                                                  */
/* Description:                                                     */
/*                                                                  */
/* Global constants I/O timer routines                              */
/*                                                                  */
/* Information:                                                     */
/*                                                                  */
/* Defines the variables required for IO timing + stats             */
/* Done globally as a number of functions will need them            */
/*                                                                  */
/* Current Code Owner: Paul Selwood                                 */
/*                                                                  */
/* History:                                                         */
/*                                                                  */
/* Version  Date       Comment                                      */
/* -------  ----       -------                                      */
/*   6.1    02/08/04   Separate deck created as part of PORTIO      */
/*                     modularisation.                  P.Dando     */
/*                                                                  */

#if defined(LINUX) || defined(DECALPHA) || defined(NEC) \
 || defined(SGI) || defined(ALTIX) || defined(IBM)
/* Linux/DEC/NEC/IBM and SGI use different includes             */
/* for gettimeofday stuff than  HP/Cray                         */
#include <sys/time.h>
#include <unistd.h>
#endif

/* Constants */
#define IO_ITER 10
#define IO_B_IN  0
#define IO_B_OUT 1
#define IO_F_OP  2
#define IO_F_CL  3

/* End C_IO_TIMER                                                   */
