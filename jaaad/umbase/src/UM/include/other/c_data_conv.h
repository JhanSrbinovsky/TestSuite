/* Start C_DATA_CONV                                                */
/*                                                                  */
/* Description:                                                     */
/*                                                                  */
/* Global constants for Data Conversion routines                    */
/*                                                                  */
/* Information:                                                     */
/*                                                                  */
/* Header file providing global constants required by the portable  */
/* data conversion routines.  These constants provide masks for     */
/* for selectively extracting the sign, exponent and mantissa       */
/* parts from IEEE and IBM 32 and 64 bit reals and integers.        */
/* The constants are expressed here in hex format.  The             */
/* corresponding binary representations are given by the usual      */
/* relationships:                                                   */
/*    0=0000, 1=0001, 2=0010,....., 9=1001, A=1010,....F=1111       */
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

/* Format 1: IBM 32 bit floating point numbers                      */
#define ibm_sign_32   0x80000000
#define ibm_expo_32   0x7F000000
#define ibm_mant_32   0x00FFFFFF
#define ibm_mant_bits_32 24
#define ibm_expo_bits_32 7
#define ibm_expo_bias_32 64

#if defined(FRL8) || defined(CRAY)
/* Format 2: IBM 64 bit floating point numbers                      */
#define ibm_sign_64   0x8000000000000000LL
#define ibm_expo_64   0x7F00000000000000LL
#define ibm_mant_64   0x00FFFFFFFFFFFFFFLL
#define ibm_mant_bits_64 56
#define ibm_expo_bits_64 7
#define ibm_expo_bias_64 64
#endif

/* Format 3: IEEE 32 bit floating point numbers                     */
#define ieee_sign_32  0x80000000
#define ieee_expo_32  0x7F800000
#define ieee_mant_32  0x007FFFFF
#define ieee_mant_bits_32 23
#define ieee_expo_bits_32 8
#define ieee_expo_bias_32 127

#if defined(FRL8) || defined(CRAY)
/* Format 4: IEEE 64 bit floating point numbers                     */
#define ieee_sign_64  0x8000000000000000LL
#define ieee_expo_64  0x7FF0000000000000LL
#define ieee_mant_64  0x000FFFFFFFFFFFFFLL
#define ieee_mant_bits_64 52
#define ieee_expo_bits_64 11
#define ieee_expo_bias_64 1023
#endif

/* Format 6: IEEE 16 bit integers (two's complement)                */
#define ieee_int_sign_16  0x8000
#define ieee_int_mant_16  0x7FFF
#define ieee_int_mant_bits_16 15

/* Format 7: IEEE 32 bit integers (two's complement)                */
#define ieee_int_sign_32  0x80000000
#define ieee_int_mant_32  0x7FFFFFFF
#define ieee_int_mant_bits_32 31

#if defined(FRL8) || defined(CRAY)
/* Format 8: IEEE 64 bit integers (two's complement)                */
#define ieee_int_sign_64  0x8000000000000000LL
#define ieee_int_mant_64  0x7FFFFFFFFFFFFFFFLL
#define ieee_int_mant_bits_64 63
#endif

/* End C_DATA_CONV                                                  */
