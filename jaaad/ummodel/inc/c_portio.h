/* Start C_PORTIO                                                   */
/*                                                                  */
/* Description:                                                     */
/*                                                                  */
/* Global definitions needed for Portable I/O                       */
/*                                                                  */
/* Information:                                                     */
/*                                                                  */
/* Provides:                                                        */
/*   1. Declarations needed to use pre-allocation via ialloc on     */
/*      Cray Research platforms.                                    */
/*   2. MAX_UNITS and BCAST global definitions                      */
/*   3. Initialises 'printstatus' control variable from its Fortran */
/*      equivalent.                                                 */
/*                                                                  */
/* Current Code Owner: Paul Dando                                   */
/*                                                                  */
/* History:                                                         */
/*                                                                  */
/* Version  Date       Comment                                      */
/* -------  ----       -------                                      */
/*   6.1    02/08/04   Separate deck created as part of PORTIO      */
/*                     modularisation.                  P.Dando     */
/*                                                                  */

#if defined(CRI_OPEN)
/* Declarations to use the pre-allocation via ialloc on Cray
   Research Platforms */
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/file.h>

#include <sys/statfs.h>
#include <sys/fstyp.h>
#include <sys/fsid.h>

#endif

#define MAX_UNITS 300

/* Define the bits used in the properies table */

#define BCAST 1

/* Initialise 'printstatus' control variable from its Fortran
   equivalent. */
#define PRINTSTATUS_MIN    1
#define PRINTSTATUS_NORMAL 2
#define PRINTSTATUS_OPER   3
#define PRINTSTATUS_DIAG   4

/* End C_PORTIO                                                     */
