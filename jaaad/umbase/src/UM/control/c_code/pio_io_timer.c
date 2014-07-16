#if defined(C95_2A) || defined(RECON) || defined(UTILIO)               \
 || defined(UTILHIST) || defined(FLDIO) || defined(VAROPSVER)
/******************************COPYRIGHT*******************************/
/* (C) Crown copyright Met Office. All rights reserved.               */
/* For further details please refer to the file COPYRIGHT.txt         */
/* which you should have received as part of this distribution.       */
/* *****************************COPYRIGHT******************************/
/**********************************************************************/
/* Description:                                                       */
/*                                                                    */
/* I/O Timing Routines                                                */
/*                                                                    */
/* Information:                                                       */
/*                                                                    */
/* Provides routines for I/O timing. The following routines are       */
/* available:                                                         */
/*                                                                    */
/*   io_eval_overhead     - Called from C.  This function attempts to */
/*                          establish how long a call to gettimeofday */
/*                          takes.                                    */
/*   io_update_timer      - Called from C. Wrapper to call the actual */
/*                          timer update function.                    */
/*   io_update_timer_data - Called from C. Function to update the     */
/*                          cumulative timers                         */
/*   io_report            - Called from C. Function to provide a      */
/*                          simple report for I/O timing.             */
/*   io_total_timings     - Called from Fortran.  Function to produce */
/*                          a simple report detailing I/O usage for   */
/*                          the complete program to date.             */
/*   io_output_inter_timings  - Called from Fortran. Function to      */
/*                              produce a simple report detailing the */
/*                              I/O usage since the last call of this */
/*                              function.                             */
/*   io_timing_init       - Called from Fortran.  Function to set the */
/*                          value of the io_timer_active flag (0 if   */
/*                          disabled; non-zero if active).            */
/*                                                                    */
/* For further details, see UMDP S5.                                  */
/*                                                                    */
/* Current Code Owner: Paul Selwood                                   */
/*                                                                    */
/* History                                                            */
/*                                                                    */
/* Version  Date       Comment                                        */
/* -------  ----       -------                                        */
/*   6.1    02/08/04   Separate deck created as part of PORTIO        */
/*                     modularisation.                      P.Dando   */
/*                                                                    */
/* Code Description:                                                  */
/*   Language: C                                                      */
/*   This code is written to UMDP3 v6 programming standards.          */
/*                                                                    */
/**********************************************************************/

/* Standard header files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Header files for Fortran to C interface */
#include "c_fort2c.h"

integer *the_unit;
static char message[256];

/* Header files for IO Timing */
#include "c_io_timer.h"

/* Global variables */
u_integer io_calls[4] = {0,0,0,0};  /* No. of calls to routine */
u_integer io_bytes[4] = {0,0,0,0};  /* Bytes transferred       */
u_integer io_secs[4]  = {0,0,0,0};  /* Seconds taken           */
u_integer io_usecs[4] = {0,0,0,0};  /* Microseconds taken      */

u_integer io_calls_last[4] = {0,0,0,0};  /* Previous values of */
u_integer io_bytes_last[4] = {0,0,0,0};  /* IO data for        */
u_integer io_secs_last[4]  = {0,0,0,0};  /* computing          */
u_integer io_usecs_last[4] = {0,0,0,0};  /* incrementals       */

u_integer io_overhead       = 0;       /* timer overhead       */
u_integer io_timer_active   = 0;       /* main IO timer switch */


/* Prototypes  for internally called functions */
void io_eval_overhead( void );
void io_update_timer( struct timeval, struct timeval,
                      u_integer, integer );
void io_update_timer_data( struct timeval, struct timeval,
                           u_integer, u_integer *,
                           u_integer *, u_integer *,
                           u_integer *);
void io_report( u_integer, u_integer, u_integer,
                u_integer, integer );


/* The following routines are for IO timing */
/****************************************************************/
/* The following function trys to make a guess at how long a    */
/* call to gettimeofday takes. It does this by making two calls */
/* immediately adjacent and taking the difference. This is      */
/* averaged over a number of calls to improve the figure given  */
/****************************************************************/
void io_eval_overhead( void )
{
  extern u_integer io_overhead;
  int i;
  struct timeval first,
                 second,
                 lapsed;

#if ! defined(LINUX)
  struct timezone tzp;
#endif

  /* if overhead is non-zero it has already been calculated. */
  if ( io_overhead != 0 )
    return;

  /* Initialised lapsed time */
  lapsed.tv_usec = 0;
  lapsed.tv_sec  = 0;

  /* Loop - calling gettimeofday the required no. of times */
  for (i=0; i<IO_ITER; i++)
  {
#if defined(LINUX)
    gettimeofday( &first,  NULL );
    gettimeofday( &second, NULL );
#else
    gettimeofday( &first,  &tzp );
    gettimeofday( &second, &tzp );
#endif

    /* Adjust for increase in no. of seconds */
    if (first.tv_usec > second.tv_usec)
    {
      second.tv_usec += 1000000;
      second.tv_sec--;
    }

    lapsed.tv_usec += ( second.tv_usec - first.tv_usec );
    lapsed.tv_sec  += ( second.tv_sec  - first.tv_sec  );

    /* Adjust for increase in no. of seconds */
    if ( lapsed.tv_usec > 1000000 )
    {
      lapsed.tv_usec -= 1000000;
      lapsed.tv_sec ++;
    }
  }

  /* Set io_overhead in usecs  */
  io_overhead = lapsed.tv_sec * 1000000;
  io_overhead += lapsed.tv_usec;

  io_overhead /= IO_ITER;

  return;
}


/****************************************************************/
/* The following function is a wrapper to call the actual timer */
/* update function for both total and last data.                */
/****************************************************************/
void io_update_timer( struct timeval first,
                      struct timeval second,
                      u_integer bytes_to_add,
                      integer op )
{
  extern u_integer io_calls[4];
  extern u_integer io_bytes[4];
  extern u_integer io_secs[4];
  extern u_integer io_usecs[4];

  io_update_timer_data( first, second, bytes_to_add,
                        &(io_calls[op]), &(io_bytes[op]),
                        &(io_secs[op]),  &(io_usecs[op]) );

  io_update_timer_data( first, second, bytes_to_add,
                        &(io_calls_last[op]), &(io_bytes_last[op]),
                        &(io_secs_last[op]),  &(io_usecs_last[op]) );
  return;
}

/****************************************************************/
/* The following function is used to update the cumulative      */
/* timers etc that are held.                                    */
/****************************************************************/
void io_update_timer_data( struct timeval first,
                           struct timeval second,
                           u_integer bytes_to_add,
                           u_integer *calls,
                           u_integer *bytes,
                           u_integer *secs,
                           u_integer *usecs)
{
  extern u_integer io_overhead;
  struct timeval lapsed;

  /* Find the lapsed time */
  if (first.tv_usec > second.tv_usec)
  {
    second.tv_usec += 1000000;
    second.tv_sec--;
  }

  lapsed.tv_usec = ( second.tv_usec - first.tv_usec );
  lapsed.tv_sec  = ( second.tv_sec  - first.tv_sec  );

  /* Update the timers */
  (*secs)  += lapsed.tv_sec;
  (*usecs) += lapsed.tv_usec;

  /* Cater for roll over of second */
  if ( *usecs  > 1000000 )
  {
    (*usecs) -= 1000000;
    (*secs) ++;
  }

  /* cater for the function call overhead */
  io_eval_overhead();

  (*usecs) -= io_overhead;

  /* increment bytes and calls */
  (*calls) ++;
  (*bytes) += bytes_to_add;

  return;
}

/****************************************************************/
/* The following function is a utility to write out a simple    */
/* report for the IO                                            */
/****************************************************************/
void io_report( u_integer calls,
                u_integer bytes,
                u_integer secs,
                u_integer usecs,
                integer op )
{
  u_integer avg_secs;
  u_integer avg_usecs;
  u_integer avg_bytes;

  integer unit;
  float time;
  float avg_time;
  float mbs;         /* Megabytes per second */
  float avg_mbs;     /* Average mbs */


  unit = -1;  /* -1 signfies unset */
  the_unit = &unit;  /* the_unit needs to be set for stdout */

  /* Calculate the averages we need for output */
  if ( calls == 0 )
  {
    avg_usecs = 0;
    avg_secs  = 0;
    avg_bytes = 0;
  }
  else
  {
    avg_usecs = ((1000000 * secs) + usecs) / calls;
    avg_secs = avg_usecs / 1000000;
    avg_usecs = avg_usecs % 1000000;

    avg_bytes = bytes/calls;
  }

  /* Convert the secs and usecs to floating point time */
  time = (float) secs + (float) usecs / 1000000.0;
  avg_time = (float) avg_secs + (float) avg_usecs / 1000000.0;

  /* do the output */
  sprintf(message, "Calls: %d", calls);
  CALL_MESSAGE_PRINT( message );

  sprintf(message,
     "Total time: %.6f secs \t\t Average time: %.6f secs",
     time, avg_time );
  CALL_MESSAGE_PRINT( message );

  if (op == IO_B_IN || op == IO_B_OUT)
  {
    sprintf(message,
           "Total bytes transferred: %ld \t Average bytes transferred:"
           " %ld", bytes, avg_bytes);
    CALL_MESSAGE_PRINT( message );

    if ( time > 0.000001 )   /* Only print transfer rates if exists */
    {
      /* Note that transfers are specified in terms of Megabytes */
      /* per second in the SI sense - ie 10^6 bytes.             */
      /* This can easily be changed for Mebibytes = 2^20 bytes   */
      mbs = (float) bytes / time;
      mbs = mbs/ (1000.0 * 1000.0);

      avg_mbs = (float) avg_bytes / avg_time;
      avg_mbs = avg_mbs/ (1000.0 * 1000.0);

      sprintf(message,
              "Overall transfer rate: %.6f MB/s \t Average transfer"
              " rate: %.6f MB/s", mbs, avg_mbs);
      CALL_MESSAGE_PRINT( message );
    }
  }

  return;
}

/****************************************************************/
/* The following function is to be called from Fortran.         */
/* It produces a simple report (on stdout) that details the     */
/* io usage for the complete program thus far.                  */
/****************************************************************/
void
#if defined(C_LOW)
io_total_timings
#elif defined(C_LOW_U)
io_total_timings_
#else
IO_TOTAL_TIMINGS
#endif
(void)
{
   extern u_integer io_calls[4];
   extern u_integer io_bytes[4];
   extern u_integer io_secs[4];
   extern u_integer io_usecs[4];

   integer unit;

   unit = -1; /* -1 is unset */
   the_unit = &unit;

   /* BUFFIN */
   sprintf(message, "\nTotal IO Timings: BUFFIN");
   CALL_MESSAGE_PRINT( message );

   io_report( io_calls[IO_B_IN], io_bytes[IO_B_IN], io_secs[IO_B_IN],
              io_usecs[IO_B_IN], IO_B_IN );


   /* BUFFOUT */
   sprintf(message, "\nTotal IO Timings: BUFFOUT");
   CALL_MESSAGE_PRINT( message );

   io_report( io_calls[IO_B_OUT], io_bytes[IO_B_OUT], io_secs[IO_B_OUT],
              io_usecs[IO_B_OUT], IO_B_OUT );


   /* FILE_OPEN */
   sprintf(message, "\nTotal IO Timings: FILE_OPEN");
   CALL_MESSAGE_PRINT( message );

   io_report( io_calls[IO_F_OP], io_bytes[IO_F_OP], io_secs[IO_F_OP],
              io_usecs[IO_F_OP], IO_F_OP );


   /* FILE_CLOSE */
   sprintf(message, "\nTotal IO Timings: FILE_CLOSE");
   CALL_MESSAGE_PRINT( message );

   io_report( io_calls[IO_F_CL], io_bytes[IO_F_CL], io_secs[IO_F_CL],
              io_usecs[IO_F_CL], IO_F_CL );

   return;
}

/****************************************************************/
/* The following function is to be called from Fortran.         */
/* It produces a simple report (on stdout) that details the     */
/* io usage since the last call of this function.               */
/****************************************************************/
void
#if defined(C_LOW)
io_output_inter_timings
#elif defined(C_LOW_U)
io_output_inter_timings_
#else
IO_OUTPUT_INTER_TIMINGS
#endif
(void)
{
   extern u_integer io_calls_last[4];
   extern u_integer io_bytes_last[4];
   extern u_integer io_secs_last[4];
   extern u_integer io_usecs_last[4];
   integer unit;
   int i;

   unit = -1; /* -1 is unset */
   the_unit = &unit;

   /* BUFFIN */
   sprintf(message, "\nIntermediate IO Timings: BUFFIN");
   CALL_MESSAGE_PRINT( message );

   io_report( io_calls_last[IO_B_IN], io_bytes_last[IO_B_IN],
              io_secs_last[IO_B_IN], io_usecs_last[IO_B_IN], IO_B_IN );


   /* BUFFOUT */
   sprintf(message, "\nIntermediate IO Timings: BUFFOUT");
   CALL_MESSAGE_PRINT( message );

   io_report( io_calls_last[IO_B_OUT], io_bytes_last[IO_B_OUT],
              io_secs_last[IO_B_OUT], io_usecs_last[IO_B_OUT],
              IO_B_OUT );


   /* FILE_OPEN */
   sprintf(message, "\nIntermediate IO Timings: FILE_OPEN");
   CALL_MESSAGE_PRINT( message );

   io_report( io_calls_last[IO_F_OP], io_bytes_last[IO_F_OP],
              io_secs_last[IO_F_OP], io_usecs_last[IO_F_OP], IO_F_OP );


   /* FILE_CLOSE */
   sprintf(message, "\nIntermediate IO Timings: FILE_CLOSE");
   CALL_MESSAGE_PRINT( message );

   io_report( io_calls_last[IO_F_CL], io_bytes_last[IO_F_CL],
              io_secs_last[IO_F_CL], io_usecs_last[IO_F_CL], IO_F_CL );

   /* reset the *last* data */
   for (i=1; i<4; i++)
   {
      io_calls_last[i] = 0;
      io_bytes_last[i] = 0;
      io_secs_last[i] = 0;
      io_usecs_last[i] = 0;
   }

   return;
}

/****************************************************************/
/* The following function is to be called from Fortran. It      */
/* is used to set the value of the io_timer_active flag.        */
/* 0        = disabled io_timer                                 */
/* non-zero = active io_timer                                   */
/****************************************************************/
void
#if defined(C_LOW)
io_timing_init
#elif defined(C_LOW_U)
io_timing_init_
#else
IO_TIMING_INIT
#endif
( integer *activate )
{
  io_timer_active = *activate;
}

#endif
