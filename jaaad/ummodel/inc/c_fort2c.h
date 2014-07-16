/* Start C_FORT2C                                                   */
/*                                                                  */
/* Description:                                                     */
/*                                                                  */
/* Definitions needed for Fortran to C Interface                    */
/*                                                                  */
/* Information:                                                     */
/*                                                                  */
/* Provides definitions needed for Fortran to C interface.          */
/* Declares variable type real to have the same wordlength          */
/* as a Fortran REAL data type, integer to have the same            */
/* wordlength as Fortran INTEGER data type.  Also provides          */
/* and interface for printing messages to stdout via the            */
/* definition CALL_MESSAGE_PRINT.                                   */
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

#if defined(FRL8) || defined(CRAY)
/* Fortran REAL is equivalent to C double */
typedef double real;
#else
/* Fortran REAL is equivalent to C float */
typedef float real;
#endif

#if defined(C_INT)
/* Fortran INTEGER is equivalent to C int */
typedef int integer;
typedef unsigned int u_integer;
#elif defined(C_LONG_INT)
/* Fortran INTEGER is equivalent to C long int */
typedef unsigned long int u_integer;
typedef long int integer;
#elif defined(C_LONG_LONG_INT)
/* Fortran INTEGER is equivalent to C long long int */
typedef unsigned long long int u_integer;
typedef long long int integer;
#else
/* DEFAULT: Fortran INTEGER is equivalent to C int */
typedef unsigned int u_integer;
typedef int integer;
#endif

/* Cray-specific Fortran characters and I/O */

#if defined(CRAY)
/* Read fortran.h to get Fortran Character stuff on Cray's */
#include <fortran.h>
#endif
#if defined(CRI_FFIO) || defined(CRI_OPEN)
/* Read Headers for FFIO */
#include <fcntl.h>
#include <ffio.h>
#endif


/* Define the function that outputs the text string for          */
/* messages.  Note that the trailing newline character is        */
/* to be supplied by the print routine, and is no longer         */
/* in the string.  Similarly, leading new lines are now handled  */
/* by the print routine - typically a newline is inserted for    */
/* each change of unit.                                          */

#if defined(CRAY)
#define CALL_MESSAGE_PRINT(text) \
 PRINT_FROM_C(the_unit, _cptofcd(text, strlen(text)))
#else
#if defined(IBM)
#if defined(C_LOW)
#define CALL_MESSAGE_PRINT(text) \
 print_from_c(the_unit,text,strlen(text))
#elif defined(C_LOW_U)
#define CALL_MESSAGE_PRINT(text) \
 print_from_c_(the_unit,text,strlen(text))
#else
#define CALL_MESSAGE_PRINT(text) \
 PRINT_FROM_C(the_unit,text,strlen(text))
#endif
#else
#define CALL_MESSAGE_PRINT(text) fprintf(stdout, "%s\n", text)
#endif
#endif

/* End C_FORT2C                                                    */
