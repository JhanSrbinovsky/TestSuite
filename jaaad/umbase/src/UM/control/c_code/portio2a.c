#if defined(C95_2A) || defined(RECON) || defined(UTILIO)               \
 || defined(UTILHIST) || defined(FLDIO) || defined(VAROPSVER)
/******************************COPYRIGHT*******************************/
/* (C) Crown copyright Met Office. All rights reserved.               */
/* For further details please refer to the file COPYRIGHT.txt         */
/* which you should have received as part of this distribution.       */
/* *****************************COPYRIGHT******************************/

/* C language routines for portable version of UM */
/*   Written by A.Dickinson 1/11/91               */
/* Standard header files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <errno.h>
#include <time.h>

/* Header files for Fortran to C interface */
#include "c_fort2c.h"
static char message[256];

/* Header files for I/O */
#include "c_portio.h"

/* Header files for I/O timing */
#include "c_io_timer.h"

/* Prototypes for internally called service functions */

extern u_integer io_timer_active ;
void io_update_timer( struct timeval, struct timeval,
                      u_integer, integer );
void change_endian( integer *, int );
void output_buffer( integer *, real *, integer *, integer *, real *);
void flush_unit_buffer( integer * );

/* Global variables used for I/O  */

#if defined(CRI_FFIO)
int fd[MAX_UNITS];
#else
FILE *pf[MAX_UNITS]=
                     {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                     };
#endif

static integer io_position[MAX_UNITS];

#if defined(BUFRD_IO)

/* Disk buffering and CRI_FFIO are not allowed together */

#if defined(CRI_FFIO)
#undef CRI_FFIO
#endif

#include <malloc.h>

/* Buffers to hold data before it is written */

static real *unit_buffer[MAX_UNITS]=
                     {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                     };

/* Pointer to the current address and offset in the buffer */

static integer unit_intent[MAX_UNITS], unit_offset[MAX_UNITS];
static integer unit_setpos[MAX_UNITS], unit_address[MAX_UNITS];

/* Define the buffer size - 32Mbytes */
#define buffer_size 4000000

/* Temporary pointer to the claimed buffer */
static char *real_add;

#endif

static integer readonly = 0;

#if defined(CRI_OPEN)

integer file_size[MAX_UNITS] =
                     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

struct statfs stats;
char buf[FSTYPSZ];

#endif
int open_flag[MAX_UNITS] =
                     {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

/* Unit properies table - one word per unit, one bit per
   property at present */

integer file_properties[MAX_UNITS] =
                     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

integer *the_unit;

/* Initialise 'printstatus' control variable.  This will
   be set from its Fortran equivalent. */
integer printstatus = 0;


#if defined(C_LOW)
extern void clear_unit_bcast_flag();
#elif defined(C_LOW_U)
extern void clear_unit_bcast_flag_();
#else
extern void CLEAR_UNIT_BCAST_FLAG();
#endif


void
#if defined(C_LOW)
set_printstatus
#elif defined(C_LOW_U)
set_printstatus_
#else
SET_PRINTSTATUS
#endif
(printstatus_f)
integer *printstatus_f;  /* IN: fortran equivalent */
{
printstatus = *printstatus_f;
}

void
#if defined(VAROPSVER)
#if defined(C_LOW)
buffin
#elif defined(C_LOW_U)
buffin_
#else
BUFFIN
#endif
#else
#if defined(C_LOW)
buffin_single
#elif defined(C_LOW_U)
buffin_single_
#else
BUFFIN_SINGLE
#endif
#endif
(unit, array, maxlen, length, status)
integer *unit;      /* Fortran unit                         */
#if defined(CRI_FFIO)
char array[];       /* Array into which data is read        */
#else
real array[];       /* Array into which data is read        */
#endif
integer *maxlen;    /* Number of real numbers to be read    */
integer *length;    /* Number of real numbers actually read */
real *status;       /* Return code                          */
{
  struct timeval first;
  struct timeval second;
#if ! defined(LINUX)
  struct timezone tzp;
#endif
#if defined(LITTLE_END)
  int i;                               /* array counter             */
  integer *ptr_integer = 0;            /* temporary integer pointer */
  void change_endian(integer *, int);  /* function to swap endian   */
#endif
integer k;
  the_unit=unit;

#if defined(CRI_FFIO)
  if(open_flag[*unit]== 0){

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
       gettimeofday( &first, &tzp );
#endif
    }

    k = ffread(fd[*unit], array, sizeof(real)*(*maxlen));

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
      gettimeofday( &second, &tzp );
#endif
      io_update_timer( first, second, k, IO_B_IN );
    }

    *length = k/sizeof(real);

        *status=-1.0;
    if(k == -1)
    {
      k=errno;
      if(k != ENOTWELLFORMED)
      {
        perror("\nBUFFIN: Read Failed");
        sprintf(message,
         "BUFFIN: C I/O Error - Return code = %d", k);
        CALL_MESSAGE_PRINT(message);
      }
      *status=5.0;
      if(k == FFEOF || k == FFEOD)
      {
        *status=0.0;
      }
      if(k == FFEOR)
      {
        *status=1.0;
      }
      if(k == FFERR)
      {
        *status=2.0;
      }
      if(k == ENOTWELLFORMED)
      {
        *status=4.0;
      }
    }
   }
   else
        *status=3.0;

    io_position[*unit]=io_position[*unit]+*length;
}

#else
  if(open_flag[*unit]== 0){

#if defined(BUFRD_IO)

/* If this file is read and write, then flush any output buffer */

    if(unit_intent[*unit] != readonly) {

      flush_unit_buffer(unit);

    }
#endif
    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
      gettimeofday( &first, &tzp );
#endif
    }

    *length = fread(array,sizeof(real),*maxlen,pf[*unit]);

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
      gettimeofday( &second, &tzp );
#endif
      io_update_timer( first, second, sizeof(real) * (*length),
                       IO_B_IN );
    }

#if defined(LITTLE_END)
    for (i = 0; i < *maxlen; i++) {
      ptr_integer = (integer *)&array[i];
      change_endian(ptr_integer, sizeof(integer));
    }
#endif

    *status=-1.0;
    k=feof(pf[*unit]);

    if(k != 0)
    {
      perror("\nBUFFIN: Read Failed");
      sprintf(message,
       "BUFFIN: C I/O Error - Return code = %d", k);
      CALL_MESSAGE_PRINT(message);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      perror("\nBUFFIN: Read Failed");
      sprintf(message,
       "BUFFIN: C I/O Error - Return code = %d", k);
      CALL_MESSAGE_PRINT(message);
      *status=1.0;
    }
   }
   else
        *status=3.0;

    io_position[*unit]=io_position[*unit]+*length;
}
#endif

void
#if defined(VAROPSVER)
#if defined(C_LOW)
buffout
#elif defined(C_LOW_U)
buffout_
#else
BUFFOUT
#endif
#else
#if defined(C_LOW)
buffout_single
#elif defined(C_LOW_U)
buffout_single_
#else
BUFFOUT_SINGLE
#endif
#endif
(unit, array, maxlen, length, status)
integer *unit;      /* Fortran unit                            */
#if defined(CRI_FFIO)
char array[];       /* Array from which data is written        */
#else
real array[];       /* Array from which data is written        */
#endif
integer *maxlen;    /* Number of real numbers to be written    */
integer *length;    /* Number of real numbers actually written */
real *status;       /* Return code                             */
{
  integer record_length;
  integer record_offset;
  integer buffer_length;
  integer i;

  integer *buffer_length_address;

  real *buffer_address;
  struct timeval first;
  struct timeval second;
#if ! defined(LINUX)
  struct timezone tzp;
#endif
  integer k;
  the_unit=unit;

#if defined(CRI_FFIO)
  if(open_flag[*unit]== 0){

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
      gettimeofday( &first, &tzp );
#endif
    }

    k = ffwrite(fd[*unit], array, sizeof(real)*(*maxlen));

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
      gettimeofday( &second, &tzp);
#endif
      io_update_timer( first, second, k, IO_B_OUT );
    }

    *length = k/sizeof(real);

        *status=-1.0;
    if(k == -1)
    {
      k=errno;
      if(k != ENOTWELLFORMED)
      {
        perror("\nBUFFOUT: Write Failed");
        sprintf(message,
         "BUFFOUT: C I/O Error - Return code = %d", k);
        CALL_MESSAGE_PRINT(message);
      }
      *status=5.0;
      if(k == FFEOF || k == FFEOD)
      {
        *status=0.0;
      }
      if(k == FFEOR)
      {
        *status=1.0;
      }
      if(k == FFERR)
      {
        *status=2.0;
      }
      if(k == ENOTWELLFORMED)
      {
        *status=4.0;
      }
    }
   }
#else

/* This is the code when CRI_FFIO is not defined */

  if(open_flag[*unit]== 0){

#if defined(BUFRD_IO)

/* Set up a few values from the input */

    record_length=*maxlen;
    record_offset=0;
    buffer_address=unit_buffer[*unit];

/* Buffered I/O - check if the new record will fit in the buffer */

    while(unit_offset[*unit]+record_length > buffer_size) {

/* Buffer is too small - store what we can and output that */

#if defined(T3E)
      i=unit_offset[*unit];
      buffer_length=buffer_size-unit_offset[*unit];

      shmem_get(&buffer_address[i], &array[record_offset],
       buffer_length, 0);

      record_offset=record_offset+buffer_length;
      record_length=record_length-buffer_length;
#else

#pragma _CRI ivdep
#pragma vdir nodep
      for(i=unit_offset[*unit]; i<buffer_size; i++) {
        buffer_address[i]=array[record_offset];
        record_length=record_length-1;
        record_offset=record_offset+1;
      }
#endif

/* Update the unit_offset */

      unit_offset[*unit]=buffer_size;
      flush_unit_buffer(unit);
    }

/* Now store the remainder of the record */

#if defined(T3E)
      i=unit_offset[*unit];
      buffer_length=record_length;

      shmem_get(&buffer_address[i], &array[record_offset],
       buffer_length, 0);

      record_offset=record_offset+buffer_length;
#else

#pragma _CRI ivdep
#pragma vdir nodep
    for(i=0; i<record_length; i++) {
      buffer_address[i+unit_offset[*unit]]=array[record_offset];
      record_offset=record_offset+1;
    }

#endif

/* Update the pointers */

    unit_offset[*unit]=unit_offset[*unit]+record_length;

    *length=*maxlen;
    *status=-1.0;
#else

    output_buffer(unit, array, maxlen, length, status);

#endif
  }
  else
       *status=3.0;

}

void output_buffer(unit, array, maxlen, length, status)
integer *unit;      /* Fortran unit                            */
#if defined(CRI_FFIO)
char array[];       /* Array from which data is written        */
#else
real array[];       /* Array from which data is written        */
#endif
integer *maxlen;    /* Number of real numbers to be written    */
integer *length;    /* Number of real numbers actually written */
real *status;       /* Return code                             */
{
integer k;
struct timeval first;
struct timeval second;
#if ! defined(LINUX)
struct timezone tzp;
#endif
integer two=2, io_time;
integer *current_position, *addr_two, *time;
#if defined(LITTLE_END)
integer i;
integer *ptr_integer;
#endif

  addr_two=&two;
  time=&io_time;

  the_unit=unit;


#if defined(LITTLE_END)
    for (i = 0; i < *maxlen; i++) {
      ptr_integer=(integer *)&array[i];
      change_endian(ptr_integer, sizeof(integer));
    }
#endif

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
      gettimeofday( &first, &tzp );
#endif
    }

    *length = fwrite(array,sizeof(real),*maxlen,pf[*unit]);

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
      gettimeofday( &second, &tzp);
#endif
      io_update_timer( first, second, sizeof(real) * (*length),
                       IO_B_OUT);
    }

#if defined(LITTLE_END)
    for (i = 0; i < *maxlen; i++) {
      ptr_integer=(integer *)&array[i];
      change_endian(ptr_integer, sizeof(integer));
    }
#endif

    *status=-1.0;
    k=feof(pf[*unit]);

    if(k != 0)
    {
      perror("\nBUFFOUT: Write Failed");
      sprintf(message,
       "BUFFOUT: C I/O Error - Return code = %d", k);
      CALL_MESSAGE_PRINT(message);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      perror("\nBUFFOUT: Write Failed");
      sprintf(message,
       "BUFFOUT: C I/O Error - Return code = %d", k);
      CALL_MESSAGE_PRINT(message);
      *status=1.0;
    }

#if defined(BUFRD_IO)

    unit_offset[*unit]=0;

#endif

    io_position[*unit]=io_position[*unit]+*length;
}

#endif
/* Test on whether we are using CRI_FFIO */


void
#if defined(VAROPSVER)
#if defined(C_LOW)
setpos
#elif defined(C_LOW_U)
setpos_
#else
SETPOS
#endif
#else
#if defined(C_LOW)
setpos_single
#elif defined(C_LOW_U)
setpos_single_
#else
SETPOS_SINGLE
#endif
#endif
(unit, word_address,err)
integer *unit;      /* Fortran unit                         */
integer *word_address; /* Number of words into file            */
integer *err;          /* Error checking, err = 0 no errors,
                                          err = 1 errors */
{
integer k;
integer byte_address;
  the_unit=unit;
#if defined(CRI_FFIO)

  if(open_flag[*unit]== 0){
    byte_address=(*word_address)*sizeof(real);
    k = ffseek(fd[*unit], byte_address, SEEK_SET);
    *err=0;
    if(k < 0){
      k=errno;
      perror("\nSETPOS: Seek Failed");
      sprintf(message,
       "SETPOS: Unit %d to Word Address %d Failed with Error Code %d",
       (int) *unit, (int) *word_address, k);
      CALL_MESSAGE_PRINT(message);
      *err=1;
      abort();
    }
  }
#else

#if defined(BUFRD_IO)
/* check if this file is for writing */

  if(unit_intent[*unit] != readonly) {

/* Does the new address match our current block address */

    if(io_position[*unit]+unit_offset[*unit] == *word_address) {

/* We already have a block position matching - no action required */

      return;

    }
/* Else we need to flush the buffer */

      flush_unit_buffer(unit);
}
#endif

    /* Only do the seek if we need to */
    if (io_position[*unit] != *word_address) {
      byte_address=(*word_address)*sizeof(real);
#if defined(LFS)
      k = fseeko(pf[*unit],byte_address,SEEK_SET);
#else /* LFS */
      k = fseek(pf[*unit],byte_address,SEEK_SET);
#endif /* LFS */
      *err = 0;
      if(k!=0){
        perror("\nSETPOS: Seek Failed");
        sprintf(message,
         "SETPOS: Unit %d to Word Address %d Failed with Error Code %d",
         (int) *unit, (int) *word_address, k);
        CALL_MESSAGE_PRINT(message);
         *err = 1;
         abort();
      }
    }
#endif
    io_position[*unit]=*word_address;

}

#if defined(BUFRD_IO)
void flush_unit_buffer(integer *unit)
{
  integer k;
  integer buffer_length;
  integer len;

  integer *buffer_length_address;
  integer *length;

  real stat;
  real minus_1=-1.0;

  real *buffer_address;
  real *status;

/* Check if there is anything to flush */

  if(unit_offset[*unit] > 0) {

    buffer_address=unit_buffer[*unit];
    buffer_length=unit_offset[*unit];
    k=buffer_length;
    buffer_length_address=&buffer_length;
    length=&len;
    status=&stat;
    output_buffer(unit, buffer_address,
                  buffer_length_address, length, status);

/* Check for an error */

    if(*status != minus_1 || *length != k) {

      fprintf(stderr,
       "\nFLUSH_UNIT_BUFFER: Error Flushing Buffered Data on PE %d\n",
       0);
      fprintf(stderr, "FLUSH_UNIT_BUFFER: Status is %7.1f\n", *status);
      fprintf(stderr,
       "FLUSH_UNIT_BUFFER: Length Requested was %d\n", k);
      fprintf(stderr,
       "FLUSH_UNIT_BUFFER: Length written   was %d\n\n", *length);
      abort();
    }
    else {
      unit_offset[*unit]=0;
    }
  }
}
#endif


void
#if defined(VAROPSVER)
#if defined(C_LOW)
file_open
#elif defined(C_LOW_U)
file_open_
#else
FILE_OPEN
#endif
#else
#if defined(C_LOW)
open_single
#elif defined(C_LOW_U)
open_single_
#else
OPEN_SINGLE
#endif
#endif
#if defined(CRAY)
(unit,f_file_name, char_len,intent,environ_var_flag,err)
_fcd f_file_name;    /* File name or environment variable    */
#else
(unit,file_name, char_len,intent,environ_var_flag,err)
char file_name[]; /* File name or environment variable    */
#endif
integer *unit;       /* Fortran unit                         */
integer *char_len;   /* No of chars in file name             */
integer *intent    ; /* =0 read only,!=0 read and write      */
integer *environ_var_flag; /* =0 file name in environment var, */
                     /*!=0 explicit file name            */
integer *err;        /* =0 file opened,                  */
                     /*!=0 file not opened               */
{
   char *fname;
#if defined(CRAY)
   char *file_name;
#endif
   struct timeval first;
   struct timeval second;
#if ! defined(LINUX)
   struct timezone tzp;
#endif
   char *gname;
#if defined(CRI_OPEN)

#define C_FLAGS O_RDWR | O_CREAT | O_RAW | O_LDRAW | O_BIG

#define I_FLAGS IA_CONT | IA_RAVL

#define J_FLAGS IA_RAVL

#define O_FLAGS O_RDWR | O_RAW | O_LDRAW | O_BIG

   struct ffsw local_stat;
   int local_fd;
   long i, k, l;
   long jj, *j;
#else
   long k;
#endif


   enum   filestat { old, new };
   enum   filestat filestatus;


   fname = calloc(*char_len + 1,1);
   the_unit=unit;
#if defined(CRI_FFIO)
   fd[*unit] = -1;
#else
   pf[*unit] = NULL;
#endif
#if defined(CRI_OPEN)
   j=&jj;
#endif
/* convert file name to C format */

#if defined(CRAY)
   file_name=__fcdtocp(f_file_name);
#endif
   strncpy( fname, file_name, *char_len );
   fname[ *char_len ] = '\0';
   sscanf( fname, "%s", fname );

   if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
     gettimeofday( &first, &tzp );
#endif
   }

   if ( *environ_var_flag == 0 )  /* File name held in environ var */
   {  gname = getenv( fname );
      if ( gname == NULL ) {
        sprintf(message,
         "OPEN:  WARNING: Environment variable %s not set",
         fname);
        CALL_MESSAGE_PRINT(message);
        open_flag[*unit]=1;
        *err=1;
        free (fname);
        return;
      }
   }
   else                           /* get file name from argmt fname */
      gname = fname;


   /* Check if file exists */

   if ( access( gname, 0 ) == 0 )  {   /* file exists */

      if (printstatus >= PRINTSTATUS_NORMAL) {
      sprintf(message,
       "OPEN:  File %s to be Opened on Unit %d Exists",
       gname, (int) *unit );
      CALL_MESSAGE_PRINT(message);
      }      /* test on printstatus */
      filestatus = old;
   }
   else  {   /* non-existent file */

      if (printstatus >= PRINTSTATUS_NORMAL) {
      sprintf(message,
       "OPEN:  File %s to be Opened on Unit %d does not Exist",
       gname, (int) *unit);
      CALL_MESSAGE_PRINT(message);
      }      /* test on printstatus */
      filestatus = new;
   }


   if ( filestatus == old )  {

      if ( *intent == readonly )  {

#if defined(CRI_FFIO)
         if(( fd[*unit] = ffopens( gname, O_RDONLY | O_RAW | O_LDRAW,
          0, 0, &local_stat, "system")) == -1 )  {
#else
         if ( ( pf[*unit] = fopen( gname, "rb" ) ) == NULL )  {
#endif

            perror("OPEN:  File Open Failed");
            sprintf(message,
              "OPEN:  Unable to Open File %s for Reading", gname );
            CALL_MESSAGE_PRINT(message);
         }
      }
      else  {   /*  *intent == read_and_write )  */

#if defined(CRI_FFIO)
         if(( fd[*unit] = ffopens( gname, O_RDWR | O_RAW | O_LDRAW,
          0, 0, &local_stat, "system")) == -1 )  {
#else
         if ( ( pf[*unit] = fopen( gname, "r+b" ) ) == NULL )  {
#endif

            perror("OPEN:  File Open Failed");
            sprintf(message,
              "OPEN:  Unable to Open File %s for Read/Write", gname );
            CALL_MESSAGE_PRINT(message);
         }
      }
   }


/* New file - check for write */
   if ( filestatus == new )  {

#if defined(CRI_FFIO)
/* Initialise the file control word */
      fd[*unit] = -1;
#else
/* Initialise the file control word to NULL */
      pf[*unit] = NULL;
#endif

      if ( *intent == readonly )  {
         sprintf(message, "OPEN:  **WARNING: FILE NOT FOUND" );
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "OPEN:  Ignored Request to Open File %s for Reading",
           gname );
         CALL_MESSAGE_PRINT(message);
      }
      else  {        /*  *intent == read_and_write   */

#if defined(CRI_OPEN)
/* Check if we should create and size the file prior to the open. */
/* This is determined by whether the size is known.               */
        if(file_size[*unit] != 0) {
/* Initially create and open the file */
          if(( local_fd = open( gname, C_FLAGS, 00755)) == -1 ) {

            perror("OPEN:  File Creation Failed");
            sprintf(message,
             "OPEN:  Unable to Open File %s for Read/Write",
             gname );
            CALL_MESSAGE_PRINT(message);
#if defined(CRI_FFIO)
            fd[*unit] = -1;
#else
            pf[*unit] = NULL;
#endif
          }
/* Successfully opened the file */
          else  {
#if !defined(VAROPSVER)
            sprintf(message, "OPEN:  File %s Created on Unit %d",
                     gname, (int) *unit );
            CALL_MESSAGE_PRINT(message);
#endif
/* Now find out about the file system we are in */
            if (fstatfs(local_fd, &stats, sizeof(struct statfs), 0)
             == -1) {
              perror("OPEN:  statfs error");
              abort();
            }

/*
            if (sysfs(GETFSTYP, stats.f_fstyp, buf) == -1) {
              perror("OPEN:  sysfs (GETFSTYP) error");
              abort();
            }

            printf("File system type = %s\n", buf);
            printf("Block size = %d\n", stats.f_bsize);
            printf("Fragment size = %d\n", stats.f_frsize);
            printf("Total number of blocks on file system = %d\n",
             stats.f_blocks);
            printf("Total number of free blocks = %d\n",
             stats.f_bfree);
            printf("Total number of file nodes (inodes) = %d\n",
             stats.f_files);
            printf("Total number of free file nodes = %d\n",
             stats.f_ffree);
            printf("Volume name = %s\n", stats.f_fname);
            printf("Pack name = %s\n", stats.f_fpack);
            printf("Primary partition bit map = %o\n",
             stats.f_priparts);
            printf("Secondary partition bit map = %o\n",
             stats.f_secparts);
            printf("Number of partitions = %d\n", stats.f_npart);
            printf("Big file threshold = %d bytes ", stats.f_bigsize);
            printf("or %d blocks\n", stats.f_bigsize/stats.f_bsize);
            printf("Big file allocation unit size = %d bytes ",
             stats.f_bigunit);
            printf("or %d blocks\n", stats.f_bigunit/stats.f_bsize);
            printf("Number of blocks in primary partitions = %d\n",
             stats.f_prinblks);
            printf("Number of free blocks in primary partitions = %d\n",
             stats.f_prinfree);
            printf("Primary partition allocation unit size = %d ",
             stats.f_priaunit);
            printf("bytes or %d blocks\n",
             stats.f_priaunit/stats.f_bsize);
            printf("Number of blocks in secondary partitions = %d\n",
             stats.f_secnblks);
            printf(
             "Number of free blocks in secondary partitions = %d\n",
             stats.f_secnfree);
            printf("Secondary partition allocation unit size = %d ",
             stats.f_secaunit);
            printf("bytes or %d blocks\n",
             stats.f_secaunit/stats.f_bsize); */

/* Round up the file size to a big file allocation unit.        */
/* This assumes the allocation unit is a correct multiple       */
/* of the allocation unit.                                      */
            i=stats.f_bigunit;
/* Check that the file allocation size is not zero - NFS files */
            if(i <= 0) i=1;
            file_size[*unit]=((file_size[*unit]+i-1)/i)*i;

/* Now extend the file to its full size */
            i=0;
            *j=0;
            l=I_FLAGS;
            k=ialloc(local_fd, file_size[*unit], l, i, j);
/* Did we managed to create a contiguous file of the right length? */
            if(k == file_size[*unit]) {

              sprintf(message,
               "OPEN:  %d Contiguous Bytes Allocated for File %s",
               (int) file_size[*unit], gname);
              CALL_MESSAGE_PRINT(message);
/* Now close the file, and reopen it as we want */
              k=close(local_fd);
              if(k != 0) {
                perror("OPEN:  Close during File Creation Failed");
                sprintf(message,
                 "OPEN:  Unable to Close Newly Created File %s",
                 gname);
                CALL_MESSAGE_PRINT(message);
                abort();
              }
/* Now open the file again */
#if defined(CRI_FFIO)
              if((fd[*unit] = ffopens(gname, O_FLAGS,
               0, 0, &local_stat, "system")) == -1)  {
#else
              if ( ( pf[*unit] = fopen( gname, "r+b" ) ) == NULL )  {
#endif
                perror("OPEN:  Re-open Newly Created File Failed");
                sprintf(message,
                 "OPEN:  Unable to Open Newly Created File %s",
                 gname);
                CALL_MESSAGE_PRINT(message);
                abort();
              }
            }
/*

  End of the Block for Contiguous Allocation */


            else {
/* Open for contiguous space failed -
   try again for non-contiguous space */
              i=0;
              *j=0;
              l=J_FLAGS;
              k=ialloc(local_fd, file_size[*unit], l, i, j);
/* Check if we have allocated space for the file */
              if(k == file_size[*unit]) {

                sprintf(message,
                 "OPEN:  %d Bytes Allocated for File %s",
                 (int) file_size[*unit], gname);
                CALL_MESSAGE_PRINT(message);
/* Now close the file, and reopen it as we want */
                k=close(local_fd);
                if(k != 0) {
                  perror("OPEN:  Close during File Creation Failed");
                  sprintf(message,
                   "OPEN:  Unable to Close Newly Created File %s",
                   gname);
                  CALL_MESSAGE_PRINT(message);
                  abort();
                }
/* Now open the file again, using fopen now */
#if defined(CRI_FFIO)
                if((fd[*unit] = ffopens(gname, O_FLAGS,
                 0, 0, &local_stat, "system")) == -1)  {
#else
                if ( ( pf[*unit] = fopen( gname, "r+b" ) ) == NULL )  {
#endif
                  perror("OPEN:  Re-open Newly Created File Failed");
                  sprintf(message,
                    "OPEN:  Unable to Open Newly Created File %s",
                   gname);
                  CALL_MESSAGE_PRINT(message);
                  abort();
                }
              }
/*

  End of the attempts to allocate space - if there seems to
  be enough space available, just create the file with no allocation. */


              else {
/* Unable to allocate space non-contiguously either - check the space
   available by finding out about the file system we are in */

                if (fstatfs(local_fd, &stats, sizeof(struct statfs), 0)
                 == -1) {
                  perror("OPEN:  statfs error");
                  abort();
                }

                if (sysfs(GETFSTYP, stats.f_fstyp, buf) == -1) {
                  perror("OPEN:  sysfs (GETFSTYP) error");
                  abort();
                }

                l=stats.f_bsize*stats.f_bfree;
/* if there is sufficient free space, use the file
   with no space allocated */
                if(l > file_size[*unit]) {
/* Now close the file, and reopen it as we want */
                  k=close(local_fd);
                  if(k != 0) {
                    perror("OPEN:  Close during File Creation Failed");
                    sprintf(message,
                     "OPEN:  Unable to Close Newly Created File %s",
                     gname);
                    CALL_MESSAGE_PRINT(message);
                    abort();
                  }
/* Now open the file again, using fopen now */
#if defined(CRI_FFIO)
                  if((fd[*unit] = ffopens(gname, O_FLAGS,
                   0, 0, &local_stat, "system")) == -1)  {
#else
                  if (( pf[*unit] = fopen( gname, "r+b" ) ) == NULL) {
#endif
                    perror("OPEN:  Re-open Newly Created File Failed");
                    sprintf(message,
                      "OPEN:  Unable to Open Newly Created File %s",
                     gname);
                    CALL_MESSAGE_PRINT(message);
                    abort();
                  }
                }
/*

  End of the attempts to pre-allocate space - there seems to be
  insufficient space available,  so we dump the stats and abort. */


                else {
/* Not enough space - pack up */
                  fprintf(stderr,
                   "\nAllocate for File %s for %d Bytes ",
                   gname, file_size[*unit]);
                  fprintf(stderr, "got Response %d\n\n", k);

                  fprintf(stderr, "File system type = %s\n", buf);
                  fprintf(stderr, "Block size = %d\n", stats.f_bsize);
                  fprintf(stderr, "Fragment size = %d\n",
                   stats.f_frsize);
                  fprintf(stderr,
                   "Total number of blocks on file system = %d\n",
                   stats.f_blocks);
                  fprintf(stderr,
                   "Total number of free blocks = %d\n",
                   stats.f_bfree);
                  fprintf(stderr,
                   "Total number of file nodes (inodes) = %d\n",
                   stats.f_files);
                  fprintf(stderr,
                   "Total number of free file nodes = %d\n",
                   stats.f_ffree);
                  fprintf(stderr, "Volume name = %s\n", stats.f_fname);
                  fprintf(stderr, "Pack name = %s\n", stats.f_fpack);
                  fprintf(stderr, "Primary partition bit map = %o\n",
                   stats.f_priparts);
                  fprintf(stderr, "Secondary partition bit map = %o\n",
                   stats.f_secparts);
                  fprintf(stderr,
                   "Number of partitions = %d\n", stats.f_npart);
                  fprintf(stderr,
                   "Big file threshold = %d bytes ", stats.f_bigsize);
                  fprintf(stderr,
                   "or %d blocks\n", stats.f_bigsize/stats.f_bsize);
                  fprintf(stderr,
                   "Big file allocation unit size = %d bytes ",
                   stats.f_bigunit);
                  fprintf(stderr,
                   "or %d blocks\n", stats.f_bigunit/stats.f_bsize);
                  fprintf(stderr,
                   "Number of blocks in primary partitions = %d\n",
                   stats.f_prinblks);
                  fprintf(stderr,
                   "Number of free blocks in primary partitions = %d\n",
                   stats.f_prinfree);
                  fprintf(stderr,
                   "Primary partition allocation unit size = %d ",
                   stats.f_priaunit);
                  fprintf(stderr, "bytes or %d blocks\n",
                   stats.f_priaunit/stats.f_bsize);
                  fprintf(stderr,
                   "Number of blocks in secondary partitions = %d\n",
                   stats.f_secnblks);
                  fprintf(stderr,
                   "No. of free blocks in secondary partitions = %d\n",
                   stats.f_secnfree);
                  fprintf(stderr,
                   "Secondary partition allocation unit size = %d ",
                   stats.f_secaunit);
                  fprintf(stderr, "bytes or %d blocks\n",
                   stats.f_secaunit/stats.f_bsize);

                  sprintf(message,
                   "OPEN:  Unable to allocate Space for File %s",
                   gname);
                  CALL_MESSAGE_PRINT(message);
                  abort();
                }
              }
            }
          }
        }
/*

  No file size given - create a zero length file */

        else {
#endif
/* File size not given - just open the file normally */
#if defined(CRI_FFIO)
          if((fd[*unit] = ffopens(gname, C_FLAGS, 00755,
           0, &local_stat, "system")) == -1)  {
#else
          if ( ( pf[*unit] = fopen( gname, "w+b" ) ) == NULL )  {
#endif

            perror("OPEN:  File Creation Failed");
            sprintf(message,
             "OPEN:  Unable to Open File %s for Read/Write", gname );
            CALL_MESSAGE_PRINT(message);
          }
          else  {
#if !defined(VAROPSVER)
            sprintf(message, "OPEN:  File %s Created on Unit %d",
                     gname, (int) *unit );
            CALL_MESSAGE_PRINT(message);
#endif
          }
#if defined(CRI_OPEN)
        }
#endif
      }
   }



   /* Set error code and open flag used by buffin and buffout */

#if defined(CRI_FFIO)
   if( fd[*unit] == -1 )  {
#else
   if( pf[*unit] == NULL )  {
#endif

      *err = 1;
      open_flag[*unit]=1;
   }
   else  {
#if defined(CRI_FFIO)
         *err = 0;
         open_flag[*unit]=0;
#else
      if ( setvbuf( pf[*unit], NULL, _IOFBF, BUFSIZ ) != 0 )  {
         perror("\n**Warning: setvbuf failed");
         *err=1;
         open_flag[*unit]=1;
       }
       else
       {
         *err = 0;
         open_flag[*unit]=0;

/*    set buffer to default size to force buffer alloc on heap */
  /*  setvbuf(pf[*unit],NULL,_IOFBF,BUFSIZ);  See above */
        }
#endif
   }
   io_position[*unit]=0;

#if defined(BUFRD_IO)

/* Check if we are just reading this unit */

  if(*intent != readonly) {

/* Check if this unit has a buffer already */

    integer malloc_error;

    if(unit_buffer[*unit] == NULL) {

/* No buffer for this file yet - create one */

/* Compute the number of bytes, including cache line alignment */

      k=(64 + buffer_size)*sizeof(real);

/* Claim memory */

      real_add=malloc(k);

/* Check the error response */

      if( real_add == NULL) {
        fprintf(stderr,
         "OPEN: PE %d had C I/O Error: failed in MALLOC\n",
         0);
        fprintf(stderr,
         "OPEN: Return code = %d, while claiming %d bytes\n",
         malloc_error, k);
        abort();
      }
      else {
        sprintf(message,
         "OPEN:  Claimed %d Bytes (%d Words) for Buffering",
         k, k/sizeof(real));
        CALL_MESSAGE_PRINT(message);
        sprintf(message,
         "OPEN:  Buffer Address is %26X", real_add);
        CALL_MESSAGE_PRINT(message);
      }

/* Adjust to a Cache Line, if we can so do */

#if defined(T3E) || defined(SV2)
      k=(long) real_add;
      k=k+64;
      k=((k >> 6) << 6);
      unit_buffer[*unit]=(real*) k;
      sprintf(message,
       "OPEN:  Buffer Address Modified to %17X for Cache",
       unit_buffer[*unit]);
      CALL_MESSAGE_PRINT(message);
#else
      unit_buffer[*unit]=(real *)real_add;
#endif

/* Set the flags for this buffer */

      unit_address[*unit]=io_position[*unit];
      unit_offset[*unit]=0;
      unit_intent[*unit]=1;
    }
  }
  else {
    unit_intent[*unit]=0;
  }
  unit_setpos[*unit]=-1;

#endif

#if defined(C_LOW)
   clear_unit_bcast_flag(unit);
#elif defined(C_LOW_U)
   clear_unit_bcast_flag_(unit);
#else
   CLEAR_UNIT_BCAST_FLAG(unit);
#endif

   if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
     gettimeofday( &second, &tzp );
#endif
     io_update_timer( first, second, 0, IO_F_OP );
   }

free (fname);
}



void
#if defined(VAROPSVER)
#if defined(C_LOW)
file_close
#elif defined(C_LOW_U)
file_close_
#else
FILE_CLOSE
#endif
#else
#if defined(C_LOW)
close_single
#elif defined(C_LOW_U)
close_single_
#else
CLOSE_SINGLE
#endif
#endif
#if defined(CRAY)
(unit,f_file_name,char_len,environ_var_flag,delete,err)
_fcd f_file_name;    /* File name or environment variable    */
#else
(unit,file_name,char_len,environ_var_flag,delete,err)
char file_name[];    /* File name or environment variable    */
#endif
integer *unit;       /* Fortran unit                         */
integer *char_len;   /* No of chars in file name             */
integer *delete;     /* =0 do not delete file,!=0 delete file*/
integer *environ_var_flag; /* =0 file name in environment var, */
                           /*!=0 explicit file name            */
integer *err; /* ERROR CHECKING err = 0 no errors, err = 1 Errors */
{
char *fname;
char *gname;
#if defined(CRAY)
char *file_name;
#endif
struct timeval first;
struct timeval second;
#if ! defined(LINUX)
struct timezone tzp;
#endif
int i;
integer k;

the_unit=unit;
fname = calloc(*char_len + 1,1);
/* first check to see if unit has been closed already (or not opened)*/
if(open_flag[*unit]== 0){    /* unit currently open  */

/* close file */
      if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
        gettimeofday( &first, &tzp );
#endif
      }

#if defined(BUFRD_IO)
      flush_unit_buffer(unit);
#endif

#if defined(CRI_FFIO)
      k=ffclose(fd[*unit]);
#else
      k=fclose(pf[*unit]);
#endif

/* convert file name to C format */
#if defined(CRAY)
        file_name=__fcdtocp(f_file_name);
#endif
        strncpy(fname,file_name,*char_len);
        fname[*char_len] = '\0';
        for (i=0; i<*char_len; i++){

            if (fname[i] == ' '){
               fname[i] = '\0';
               break;
            }
         }

        if(*environ_var_flag == 0)
          { gname = getenv( fname );
            if ( gname == NULL ) {
              sprintf(message,
               "CLOSE: WARNING: Environment variable %s not set",
               fname);
              CALL_MESSAGE_PRINT(message);
            open_flag[*unit]=1;
            *err=1;
            free( fname );
            return;
            }
          }
        else
          gname=fname;

#if defined(CRI_FFIO)
       if(k >= 0){
#else
      if(k==0){
#endif

/* delete file */
        if(*delete != 0){
          k=remove(gname);
          if( k != 0){
            sprintf(message,
             "CLOSE: Cannot Delete File %s",gname);
            CALL_MESSAGE_PRINT(message);
            *err = 1;
            abort();
          }
          else{  /*normal end to delete so:*/
            open_flag[*unit]=1;     /* set unit flag to closed */
            *err = 0;
            sprintf(message,
             "CLOSE: File %s Deleted", gname);
            CALL_MESSAGE_PRINT(message);
          }

        }
        else{
/* file closed */
           open_flag[*unit]=1;     /* set unit flag to closed */
          sprintf(message,
           "CLOSE: File %s Closed on Unit %d",
           gname, (int) *unit);
          CALL_MESSAGE_PRINT(message);
        }
      }
/* file not closed */
    else {
          sprintf(message,
           "CLOSE: Cannot Close File %s on Unit %d",
           gname, (int) *unit);
          CALL_MESSAGE_PRINT(message);
    }

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
      gettimeofday( &second, &tzp );
#endif
      io_update_timer( first, second, 0, IO_F_CL );
    }

}   /* end of test for unit already closed */

else {
/* unit either closed already or not open yet */
      if (printstatus >= PRINTSTATUS_NORMAL) {
  sprintf(message,
   "CLOSE: WARNING: Unit %d Not Opened", (int) *unit);
  CALL_MESSAGE_PRINT(message);
      }      /* test on printstatus */
}

free( fname );
}


void
#if defined(C_LOW)
date_time
#elif defined(C_LOW_U)
date_time_
#else
DATE_TIME
#endif
(year,month,day,hour,minute,second)
integer *year;       /* year                  */
integer *month;      /* month   Current date  */
integer *day;        /* day      and time     */
integer *hour;       /* hour                  */
integer *minute;     /* minute                */
integer *second;     /* second                */

{
char s[5];
time_t t,*r,a;
    r=&a;
    t=time(r);
    strftime(s,5,"%Y",localtime(r));
    *year=atoi(s);
    strftime(s,5,"%m",localtime(r));
    *month=atoi(s);
    strftime(s,5,"%d",localtime(r));
    *day=atoi(s);
    strftime(s,5,"%H",localtime(r));
    *hour=atoi(s);
    strftime(s,5,"%M",localtime(r));
    *minute=atoi(s);
    strftime(s,5,"%S",localtime(r));
    *second=atoi(s);
}
void
#if defined(C_LOW)
shell
#elif defined(C_LOW_U)
shell_
#else
SHELL
#endif
#if defined(CRAY)
(f_command,command_len)
_fcd f_command  ;    /* Command to be executed               */
#else
(command,command_len)
char command  []; /* Command to be executed               */
#endif
integer *command_len; /* No of chars in command               */
{
char *fname;
int i;
#if defined(CRAY)
char *command;
#endif

/* convert file name to C format */

fname = calloc( *command_len + 1,1);
#if defined(CRAY)
command=__fcdtocp(f_command);
#endif
strncpy(fname,command,*command_len);
fname[*command_len]='\0';

/* execute command */
        i=system(fname);
        if ( i == -1 ){
          /* command failed */
          printf("C Error: failed in SHELL\n");
        }
free( fname );
}

void
#if defined(VAROPSVER)
#if defined(C_LOW)
buffo32
#elif defined(C_LOW_U)
buffo32_
#else
BUFFO32
#endif
#else
#if defined(C_LOW)
buffo32_single
#elif defined(C_LOW_U)
buffo32_single_
#else
BUFFO32_SINGLE
#endif
#endif
(unit, array, maxlen, length, status)
integer *unit;     /* Fortran unit                            */
#if defined(CRI_FFIO)
char array[];      /* Array from which data is written        */
#else
real array[];      /* Array from which data is written        */
#endif
integer *maxlen;   /* Number of real numbers to be written    */
integer *length;   /* Number of real numbers actually written */
real *status;      /* Return code                             */
{
  struct timeval first;
  struct timeval second;
#if ! defined(LINUX)
  struct timezone tzp;
#endif
#if defined(LITTLE_END)
    int i;                              /* array counter             */
    integer *ptr_integer = 0;           /* temporary integer pointer */
    void change_endian(integer *,int);  /* function to swap endian   */
#endif
integer k;

    if (open_flag[*unit]==0){
#if defined(CRI_FFIO)

    if (io_timer_active) {
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
      gettimeofday( &first, &tzp );
#endif
    }

    k = ffwrite(fd[*unit], array, 4*(*maxlen));

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
      gettimeofday( &second, &tzp);
#endif
      io_update_timer( first, second, k, IO_B_OUT );
    }

    *length = k/4;

        *status=-1.0;
    if(k == -1)
    {
      k=errno;
      printf("C I/O Error: failed in BUFFO32\n");
      printf("Return code = %d\n",k);
      if(k == FFEOF || k == FFEOD)
      {
        *status=0.0;
      }
      if(k == FFEOR)
      {
        *status=1.0;
      }
      if(k == FFERR)
      {
        *status=2.0;
      }
#else

#if defined(LITTLE_END)
#if defined (FRL8) || defined (CRAY) || defined (IEEE)
     for (i = 0; i < (*maxlen + 1)/2; i++) {
#else
     for (i = 0; i < *maxlen; i++) {
#endif
        ptr_integer = (integer *)&array[i] ;
#if defined (IEEE)
        change_endian(ptr_integer, sizeof(integer));
#else
        change_endian(ptr_integer, 4);
#endif
      }
#endif

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
      gettimeofday( &first, &tzp );
#endif
    }

#if defined(BUFRD_IO)
/* As this subroutine does not perform buffered IO, but     */
/* may use an output unit that has already been buffered,   */
/* we need to flush the existing buffer before continuing.  */

    flush_unit_buffer(unit);
#endif

    *length = fwrite(array,4,*maxlen,pf [*unit]);

    if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
      gettimeofday( &second, &tzp);
#endif
      io_update_timer( first, second, 4 * (*maxlen),
                       IO_B_OUT);
    }

#if defined(LITTLE_END)
#if defined (FRL8) || defined (CRAY) || defined (IEEE)
      for (i = 0; i < (*maxlen + 1)/2; i++) {
#else
      for (i = 0; i < *maxlen; i++) {
#endif
        ptr_integer = (integer *)&array[i] ;
#if defined (IEEE)
        change_endian(ptr_integer, sizeof(integer));
#else
        change_endian(ptr_integer, 4);
#endif
      }
#endif

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFO32\n");
      printf("Return code = %d\n",k);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFO32\n");
      printf("Return code = %d\n",k);
      *status=1.0;
#endif
      }
    }
      else
        *status=3.0;

#if defined(BUFRD_IO)
/* the position in the file must be updated */
if (*length % 2 != 0)
{
  fprintf(stderr,
 "\nWARNING BUFFO32_SINGLE and BUFRD_IO not compatible\n");
}
io_position[*unit]=io_position[*unit]+*length/2;
#endif


}

void
#if defined(VAROPSVER)
#if defined(C_LOW)
buffin32
#elif defined(C_LOW_U)
buffin32_
#else
BUFFIN32
#endif
#else
#if defined(C_LOW)
buffin32_single
#elif defined(C_LOW_U)
buffin32_single_
#else
BUFFIN32_SINGLE
#endif
#endif
(unit, array, maxlen, length, status)
integer *unit;     /* Fortran unit                         */
#if defined(CRI_FFIO)
char array[];      /* Array from which data is read        */
#else
real array[];     /* Array from which data is read        */
#endif
integer *maxlen;   /* Number of real numbers to be read    */
integer *length;   /* Number of real numbers actually read */
real *status;      /* Return code                          */
{
  struct timeval first;
  struct timeval second;
#if ! defined(LINUX)
  struct timezone tzp;
#endif
#if defined(LITTLE_END)
    int i;                              /* array counter             */
    integer *ptr_integer = 0;           /* temporary integer pointer */
    void change_endian(integer *,int);  /* function to swap endian   */
#endif
integer k;

#if defined(BUFRD_IO)
/* As this subroutine does not perform buffered IO, but      */
/* may use a unit that has already been buffered, we need to */
/* flush the existing buffer before continuing to ensure we  */
/* read up to date data.                                     */

    flush_unit_buffer(unit);
#endif

    if (open_flag[*unit]==0){

#if defined(CRI_FFIO)
      if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
         gettimeofday( &first, &tzp );
#endif
      }

      k = ffread(fd[*unit], array, 4*(*maxlen));

      if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
        gettimeofday( &second, &tzp );
#endif
        io_update_timer( first, second, k, IO_B_IN );
      }

      *length = k/4;

      *status=-1.0;
      if(k == -1)
      {
        k=errno;


        if(k == FFEOF || k == FFEOD)
        {
          *status=0.0;
        }
        if(k == FFEOR)
        {
          *status=1.0;
        }
        if(k == FFERR)
        {
          *status=2.0;
        }
#else
      if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &first, NULL );
#else
         gettimeofday( &first, &tzp );
#endif
      }

      *length = fread(array,4,*maxlen,pf [*unit]);

      if (io_timer_active){
#if defined(LINUX)
       gettimeofday( &second, NULL );
#else
        gettimeofday( &second, &tzp );
#endif
        io_update_timer( first, second, 4 * (*maxlen),
                         IO_B_IN );
      }

#if defined(LITTLE_END)
#if defined (FRL8) || defined (CRAY)
      for (i = 0; i < (*maxlen + 1)/2; i++) {
#else
      for (i = 0; i < *maxlen; i++) {
#endif
        ptr_integer = (integer *)&array[i];
        change_endian(ptr_integer, 4);
      }
#endif

      *status=-1.0;
      k=feof(pf[*unit]);
      if(k != 0)
      {
        printf("C I/O Error: failed in BUFFIN32\n");
        printf("Return code = %d\n",k);
        *status=0.0;
      }
      k=ferror(pf[*unit]);
      if(k != 0)
      {
        printf("C I/O Error: failed in BUFFIN32\n");
        printf("Return code = %d\n",k);
        *status=1.0;
#endif
      }
    }
    else
      *status=3.0;
#if defined(BUFRD_IO)
/* the position in the file must be updated */
    if (*length % 2 != 0)
    {
      fprintf(stderr,
     "\nWARNING BUFFIN32_SINGLE and BUFRD_IO not compatible\n");
    }
    io_position[*unit]=io_position[*unit]+*length/2;
#endif
}

void
#if defined(VAROPSVER)
#if defined(C_LOW)
buffin8
#elif defined(C_LOW_U)
buffin8_
#else
BUFFIN8
#endif
#else
#if defined(C_LOW)
buffin8_single
#elif defined(C_LOW_U)
buffin8_single_
#else
BUFFIN8_SINGLE
#endif
#endif
#if defined(CRAY)
(unit, f_array, maxlen, length, status)
#else
(unit, array, maxlen, length, status)
#endif
integer *unit;     /* Fortran unit                         */
#if defined(CRAY)
_fcd f_array; /* Array into which data is read        */
#else
char array[];      /* Array into which data is read        */
#endif
integer *maxlen;   /* Number of bytes to be read           */
integer *length;   /* Number of bytes actually read        */
real *status;      /* Return code                          */
{
integer k;
#if defined(CRAY)
char *array;

  array=__fcdtocp(f_array);
#endif

#if defined(CRI_FFIO)
  if(open_flag[*unit]== 0){
    *length = ffread(fd[*unit], array, *maxlen);

        *status=-1.0;
    if(*length == -1)
    {
      k=errno;
      printf("C I/O Error: failed in BUFFIN8\n");
      printf("Return code = %d\n",k);
      if(k == FFEOF || k == FFEOD)
      {
        *status=0.0;
      }
      if(k == FFEOR)
      {
        *status=1.0;
      }
      if(k == FFERR)
      {
        *status=2.0;
      }
    }
   }
#else
  if(open_flag[*unit]== 0){
    *length = fread(array,1,*maxlen,pf[*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFIN8\n");
      printf("Return code = %d\n",k);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFIN8\n");
      printf("Return code = %d\n",k);
      *status=1.0;
    }
   }
#endif
   else
        *status=3.0;

}

void
#if defined(C_LOW)
buffou8
#elif defined(C_LOW_U)
buffou8_
#else
BUFFOU8
#endif
#if defined(CRAY)
(unit, f_array, maxlen, length, status)
#else
(unit, array, maxlen, length, status)
#endif
integer *unit;     /* Fortran unit                            */
#if defined(CRAY)
_fcd f_array; /* Array into which data is read        */
#else
char   array[];    /* Array from which data is written        */
#endif
integer *maxlen;   /* Number of bytes to be written           */
integer *length;   /* Number of bytes actually written        */
real *status;      /* Return code                             */
{
integer k;
#if defined(CRAY)
char *array;

  array=__fcdtocp(f_array);
#endif

#if defined(CRI_FFIO)
  if(open_flag[*unit]== 0){
    *length = ffwrite(fd[*unit], array, *maxlen);

        *status=-1.0;
    if(*length == -1)
    {
      k=errno;
      printf("C I/O Error: failed in BUFFOU8\n");
      printf("Return code = %d\n",k);
      if(k == FFEOF || k == FFEOD)
      {
        *status=0.0;
      }
      if(k == FFEOR)
      {
        *status=1.0;
      }
      if(k == FFERR)
      {
        *status=2.0;
      }
    }
   }
#else
  if(open_flag[*unit]== 0){
    *length = fwrite(array,1,*maxlen,pf[*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFOU8\n");
      printf("Return code = %d\n",k);
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      printf("C I/O Error: failed in BUFFOU8\n");
      printf("Return code = %d\n",k);
      *status=1.0;
    }
   }
#endif
   else
        *status=3.0;

}
void
#if defined(VAROPSVER)
#if defined(C_LOW)
setpos8
#elif defined(C_LOW_U)
setpos8_
#else
SETPOS8
#endif
#else
#if defined(C_LOW)
setpos8_single
#elif defined(C_LOW_U)
setpos8_single_
#else
SETPOS8_SINGLE
#endif
#endif
(unit, byte_address)
integer *unit;     /* Fortran unit                         */
integer *byte_address; /* Number of bytes into file            */
{
integer k;

#if defined(CRI_FFIO)
    k = ffseek(fd[*unit], *byte_address, SEEK_SET);
    if(k < 0){
      k=errno;
      printf("ERROR detected in SETPOS\n");
      printf("word_address = %d\n", (int) *byte_address);
      printf("Return code from fseek = %d\n",k);
    }
#else
#if defined(LFS)
    k = fseeko(pf[*unit],*byte_address,SEEK_SET);
#else /* LFS */
    k = fseek(pf[*unit],*byte_address,SEEK_SET);
#endif /* LFS */
    if(k!=0){
         printf("ERROR detected in SETPOS\n");
         printf("word_address = %d\n", (int) *byte_address);
         printf("Return code from fseek = %d\n",k);
    }
#endif
}

void
#if defined(C_LOW)
getpos8
#elif defined(C_LOW_U)
getpos8_
#else
GETPOS8
#endif
(unit, byte_address)
integer *unit;     /* Fortran unit                         */
integer *byte_address; /* Number of bytes into file            */
{

#if defined(CRI_FFIO)
    printf("Illegal Call to GETPOS8 for Unit %d\n",*unit);
    abort();
#else
#if defined(LFS)
    *byte_address = ftello(pf[*unit]);
#else /* LFS */
    *byte_address = ftell(pf[*unit]);
#endif /* LFS */
#endif


}
void
#if defined(C_LOW)
setpos32
#elif defined(C_LOW_U)
setpos32_
#else
SETPOS32
#endif
(unit,word32_address,err)

integer *unit;            /* Fortran unit                         */
integer *word32_address;  /* Number of 32bit words into the file  */
integer *err;             /* 0: no error
                             1: error occured                     */

{
integer k;
integer byte_address;

  *err=0;
  if (open_flag[*unit]==0) {
    byte_address=(*word32_address)*4;

#if defined(CRI_FFIO)
    k=ffseek(fd[*unit],byte_address,SEEK_SET);
#else
#if defined(LFS)
    k=fseeko(pf[*unit],byte_address,SEEK_SET);
#else /* LFS */
    k=fseek(pf[*unit],byte_address,SEEK_SET);
#endif /* LFS */
#endif

    if (k != 0) {
      k=errno;
      perror("\nSETPOS32: Seek failed");
      sprintf(message,
        "SETPOS32: Unit %d to 32bit Word Address %d failed. Error: %d",
        (int) *unit, (int) *word32_address, k);
      the_unit=unit;
      CALL_MESSAGE_PRINT(message);
      *err=1;
      abort();
    }

  }
}

void
#if defined(C_LOW)
getpos32
#elif defined(C_LOW_U)
getpos32_
#else
GETPOS32
#endif
(unit,word32_address)

integer *unit;            /* Fortran unit                         */
integer *word32_address;  /* Number of 32bit words into the file  */

{
int byte_address;

  if (open_flag[*unit]==0) {
#if defined(CRI_FFIO)
    printf("Illegal Call to GETPOS8 for Unit %d\n",*unit);
    abort();
#else
#if defined(LFS)
    byte_address=ftello(pf[*unit]);
#else /* LFS */
    byte_address=ftell(pf[*unit]);
#endif /* LFS */
    *word32_address = byte_address/4;
#endif
  }
}

void
#if defined(C_LOW)
getpos
#elif defined(C_LOW_U)
getpos_
#else
GETPOS
#endif
(unit, word_address)
integer *unit;     /* Fortran unit                         */
integer *word_address; /* Number of words into file            */
{
long byte_address;

   the_unit=unit;
#if defined(CRI_FFIO)
    *word_address = io_position[*unit];
/*   sprintf(message,
    "GETPOS: Illegal Call to GETPOS for Unit %d", (int) *unit);
   CALL_MESSAGE_PRINT(message); */
#else
#if defined(LFS)
     byte_address = ftello(pf[*unit]);
#else /* LFS */
     byte_address = ftell(pf[*unit]);
#endif /* LFS */
    *word_address = byte_address/sizeof(real);
#endif
      if(*word_address != io_position[*unit]) {
        sprintf(message,
         "GETPOS: IO_POSITION is %d, but FTELL gives %d",
         (int) io_position[*unit], (int) *word_address);
        the_unit=unit;
        CALL_MESSAGE_PRINT(message);
        abort();
      }

}

void
#if defined(C_LOW)
word_length
#elif defined(C_LOW_U)
word_length_
#else
WORD_LENGTH
#endif
(length)
integer *length;  /* Word length used by hardware         */
{

    *length=sizeof(real);

}
void
#if defined(C_LOW)
get_file
#elif defined(C_LOW_U)
get_file_
#else
GET_FILE
#endif
#if defined(CRAY)
(unit,f_filename,file_len,err)
_fcd f_filename; /* File name                           */
#else
(unit,filename,file_len,err)
char filename[]; /* File name                           */
#endif
integer *file_len; /* Dimension of filename               */
integer  *unit;    /* Fortran unit number                 */
integer  *err; /* Error checking err = 0 no errors, err = 1 errors */
{
char fname[16];
char fno[4];
char *gname;
int i;
int k;
#if defined(CRAY)
char *filename;
#endif

the_unit=unit;
#if defined(CRAY)
filename=__fcdtocp(f_filename);
#endif
/* construct environment variable name         */
/* in form UNITnn, where nn is Fortran unit no */

       if ( *unit < 100){
       strcpy (fname, "UNIT");
       sprintf(fno,"%02i", (int) *unit);
       strcat(fname,fno);
       fname[6]='\0';
       }
       else{
       strcpy (fname, "UNIT");
       sprintf(fno,"%03i", (int) *unit);
       strcat(fname,fno);
       fname[7]='\0';
       }

/* get file name stored in environment variable UNITnn */
       gname=getenv(fname);
       if ( gname == NULL) {
         sprintf(message,
          "GET_FILE: Environment Variable %s not Set", fname);
         CALL_MESSAGE_PRINT(message);
         filename[0] = '\0';
         for (i=1; i<*file_len; i++){
                filename[i] = ' ';
         }
         return;
       }
       k=strlen(gname);
       if(k >  *file_len){
         sprintf(message,
          "GET_FILE: File Name too long for Allocated Storage");
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "GET_FILE: Environment Variable %s", fname);
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "GET_FILE: File Name %s", gname);
         CALL_MESSAGE_PRINT(message);
         *err = 1;
         abort();
       }

/* convert file name to Fortran format */
          *err = 0;
          strncpy(filename,gname,k);
          for (i=k; i<*file_len; i++){
               filename[i] = ' ';
           }


}
void
#if defined(C_LOW)
abort()
#elif defined(C_LOW_U)
abort_()
#else
#if defined(CRAY_MSG)
ABORT(_fcd f_msg)
#else
ABORT()
#endif
#endif
{
#if defined(CRAY_MSG)
  char *msg=__fcdtocp(f_msg);
  int msg_len=__fcdlen(f_msg);
  char ch_msg[msg_len+1];
  integer *unit, *iostat;
  integer six=6, zero;

  unit=&six;
  iostat=&zero;
  FLUSH(unit, iostat);
  if(*iostat != 0)
  {
  fprintf(stderr,
   "\nUM_ABORT: Return Value from FLUSH was %d\n", (int) *iostat);
  }

  strncpy(ch_msg, msg, msg_len);
  ch_msg[msg_len]='\0';
  printf("\nUM ABORT: %s\n", ch_msg);
#endif
  abort();
}
     /*          Force i/o buffer to be written to file */
     /*          explicitly to prevent continuation run */
     /*          problems following 'hard' failures.    */
void
#if defined(C_LOW)
flush_buffer
#elif defined(C_LOW_U)
flush_buffer_
#else
FLUSH_BUFFER
#endif
(unit, icode)
integer  *unit     ;  /* Fortran unit number             */
integer  *icode    ;  /* Integer return code             */
{
int  i         ;
#if defined(CRI_FFIO)
  if(open_flag[*unit]== 0){
      i =   ffflush(fd[*unit]);
      *icode = i;
  }
#else
  if(open_flag[*unit]== 0){

#if defined(BUFRD_IO)
      flush_unit_buffer(unit);
#endif


      i =   fflush(pf[*unit]);
      *icode = i;
  }
  else {
    if(pf[*unit] != NULL) {
      sprintf(message,
       "FLUSH_BUFFER: File Pointer for Unopened Unit %d is %16X",
       (int) *unit, (unsigned long) pf[*unit]);
      the_unit=unit;
      CALL_MESSAGE_PRINT(message);
      abort();
    }
  }
#endif
}
void
#if defined(C_LOW)
fort_get_env
#elif defined(C_LOW_U)
fort_get_env_
#else
FORT_GET_ENV
#endif
#if defined(CRAY)
(f_env_var_name,ev_len,f_ev_contents,cont_len,ret_code)
_fcd f_env_var_name; /* Name of environment variable   */
integer *ev_len;     /* length of name                 */
_fcd f_ev_contents;  /* contents of environment variable */
#else
(env_var_name,ev_len,ev_contents,cont_len,ret_code)
char *env_var_name;  /* Name of environment variable   */
integer *ev_len;     /* length of name                 */
char *ev_contents;   /* contents of environment variable */
#endif
integer *cont_len;   /* length of contents              */
integer *ret_code;   /* return code: 0=OK  -1=problems  */

{
        integer minus_one=-1;

#if defined(CRAY)
        char *env_var_name=__fcdtocp(f_env_var_name);
        char *ev_contents=__fcdtocp(f_ev_contents);
#endif
        char *c_env_var_name;

        char *value;
        int len,i;

        c_env_var_name = calloc(*ev_len + 1,1);
        the_unit=&minus_one;
        strncpy(c_env_var_name,env_var_name,*ev_len);
        c_env_var_name[*ev_len]='\0';
        sscanf(c_env_var_name,"%s",c_env_var_name);

        value=getenv(c_env_var_name);
        if (value==NULL){
          *ret_code=-1;
          free( c_env_var_name );
          return;}
        else{
          *ret_code=0;}

        len=strlen(value);
        if (len > *cont_len){
         sprintf(message,
          "FORT_GET_ENV: Value too long for Allocated Storage");
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "FORT_GET_ENV: Environment Variable %s",
          c_env_var_name);
         CALL_MESSAGE_PRINT(message);
         sprintf(message,
          "FORT_GET_ENV: Value %s", value);
         CALL_MESSAGE_PRINT(message);
                abort();
        }

        strncpy(ev_contents,value,len);
        for (i=len; i<*cont_len; i++){
                ev_contents[i]=' ';
        }
free( c_env_var_name );

}
#if defined(CRI_FFIO)

void
CLOSE_ALL_FILES
()
{
integer minus_one=-1;
  int unit, k;

    sprintf(message,
     "\nCall to Close All Files");
    CALL_MESSAGE_PRINT(message);
    sprintf(message,
     "-----------------------\n");
    CALL_MESSAGE_PRINT(message);
/* Loop over all known Units */
    for (unit=0; unit<MAX_UNITS; unit++)
    {
/* first check to see if unit has been closed already (or not opened)*/
        if(open_flag[unit]== 0)
        {
            the_unit=&minus_one;
/* close file */
            k=ffclose(fd[unit]);
/* check the error response */
            if(k >= 0)
            {
/* set unit flag to closed */
                open_flag[unit]=1;
                sprintf(message,
                 "CLOSE: File on Unit %3d Closed", (int) unit);
                CALL_MESSAGE_PRINT(message);
            }
/* file not closed */
            else
            {
                sprintf(message,
                 "CLOSE: Cannot Close File on Unit %3d - Code = %d",
                 (int) unit, errno);
                CALL_MESSAGE_PRINT(message);
            }
        }
    }
    sprintf(message, "\n");
    CALL_MESSAGE_PRINT(message);
}
void
FLUSH_ALL_FILES
()
{
integer minus_one=-1;
  int unit, k;

    sprintf(message,
     "\nCall to Flush All Files");
    CALL_MESSAGE_PRINT(message);
    sprintf(message,
     "-----------------------\n");
    CALL_MESSAGE_PRINT(message);
/* Loop over all known Units */
    for (unit=0; unit<MAX_UNITS; unit++)
    {
/* first check to see if unit has been Flushed already (or not opened)*/
        if(open_flag[unit]== 0)
        {
            the_unit=&minus_one;
/* Flush file */
            k=ffflush(fd[unit]);
/* check the error response */
            if(k >= 0)
            {
/* set unit flag to Flushed */
                sprintf(message,
                 "FLUSH: File on Unit %3d Flushed", (int) unit);
                CALL_MESSAGE_PRINT(message);
            }
/* file not flushed */
            else
            {
                sprintf(message,
                 "FLUSH: Cannot Flush File on Unit %3d - Code = %d",
                 (int) unit, errno);
                CALL_MESSAGE_PRINT(message);
            }
        }
    }
    sprintf(message, "\n");
    CALL_MESSAGE_PRINT(message);
}
#else

void
CLOSE_ALL_FILES
()
{
}

void
FLUSH_ALL_FILES
()
{
}

#endif
/*                                                              */
/* Entry to accept the current File length prior                */
/* to an open request                                           */
/*                                                              */
void
#if defined(C_LOW)
set_dumpfile_length
#elif defined(C_LOW_U)
set_dumpfile_length_
#else
SET_DUMPFILE_LENGTH
#endif
(unit, length)
integer *unit;
integer *length;
{
#if defined(CRI_OPEN)
  file_size[*unit]=*length*sizeof(real);
  the_unit=unit;
  sprintf(message,
   "File Length for Unit %d set to %d Bytes",
   (int) *unit, (int) file_size[*unit]);
  CALL_MESSAGE_PRINT(message);
#endif
}
void
#if defined(C_LOW)
clear_unit_bcast_flag
#elif defined(C_LOW_U)
clear_unit_bcast_flag_
#else
CLEAR_UNIT_BCAST_FLAG
#endif
(unit)
integer *unit;
{
  integer minus_one=-1;

  file_properties[*unit]=(minus_one ^ BCAST) & file_properties[*unit];

}

void
#if defined(C_LOW)
set_unit_bcast_flag
#elif defined(C_LOW_U)
set_unit_bcast_flag_
#else
SET_UNIT_BCAST_FLAG
#endif
(unit)
integer *unit;
{

  file_properties[*unit]=BCAST | file_properties[*unit];

}


void
#if defined(C_LOW)
find_unit_bcast_flag
#elif defined(C_LOW_U)
find_unit_bcast_flag_
#else
FIND_UNIT_BCAST_FLAG
#endif
(unit, flag)
integer *unit;
integer *flag; /* non-zero if set, otherwise 0 */
{

  *flag=BCAST & file_properties[*unit];

}


#if defined(LITTLE_END)

/* In order to be consistent with present systems all data files      */
/* will be saved in big endian format.                                */
/*                                                                    */
/* The following functions are used to convert byte order (endian)    */
/* of data                                                            */
/*                                                                    */
/* Use only when reading (writing) data on little endian machines     */

void change_endian ( integer *ptr_array, int Nbytes ) {

/* Swap byte order of *ptr_array */

  int i;
  unsigned char *ptr_IVal=0;
  unsigned char *ptr_OVal=0;

  integer IVal=*(integer *)ptr_array ;
  ptr_IVal=(unsigned char *)&IVal ;

  /* unsigned char is one byte */

  ptr_OVal=(unsigned char *)ptr_array;

  for (i=0; i<Nbytes; i++) {
    /* reverse byte ordering */
    ptr_OVal[Nbytes-1-i]=ptr_IVal[i];
  }

  /* Packed data (Nbytes=4) needs the second 4 bytes swapping too */

  if ( Nbytes == 4 && Nbytes != sizeof(integer) ) {
    for (i=0; i<Nbytes; i++){
      /* reverse byte ordering */
      ptr_OVal[2*Nbytes-1-i]=ptr_IVal[Nbytes+i];
    }
  }
  return;
}

#endif

/****************************************************************/
/* The following function is to be called from Fortran.         */
/* It determines the status (open or closed) of a unit.         */
/* Returned values are C style ie,                              */
/*                     0    -    unit is open                   */
/*                     1    -    unit is closed                 */
/****************************************************************/
void
#if defined(C_LOW)
is_unit_open
#elif defined(C_LOW_U)
is_unit_open_
#else
IS_UNIT_OPEN
#endif
(unit, ret_code)
integer *unit;
integer *ret_code; /* 0 if open - 1 if closed (C style) */
{

  *ret_code = open_flag[*unit];
  return;
}
#endif
