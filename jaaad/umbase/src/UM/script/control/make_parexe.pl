#!/usr/bin/perl
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
#  Script: make_parexe.pl
#
#  Purpose:  Create a parexe file which provides all the environment
#            as well as the module to run. This is then "mpirun"ed
#            to ensure all CPUs have the full envrionment available 
#            to them.
#
#            Implementation in perl to ensure IO is done in a single
#            block and is thus efficient.
#
#  Current Owner: P Selwood      Date: 25/09/06
#  Reviewer: M Saunby
#
#   Tested under perl v5.6.1: unknown
#   UM version no: 6.1 and 6.3
#
#   History:
#    Model
#  version  Date         Modification history:
#  6.3      25/09/06     Initial version. P.Selwood.
#  6.4      07/02/07     To run on Linux PC, ulimit needs setting
#                        S.D.Mullerworth
#  6.4      21/02/07     Reversing previous change to set ulimit.
#
# External documentation: none
#
# Interface and arguments:
#
#  make_parexe.pl <none>
#
# Imports:
#   PAREXE
#   PWD
#   LOADMODULE
#
# Exports: none
#
# End of header -------------------------------------------------------

# find the PAREXE environment variable and open the file referred to.
$file = $ENV{"PAREXE"} 
      or die "make_parexe: Couldn't evaluate variable PAREXE\n";
open (PF, ">$file") 
      or die "make_parexe: Couldn't open file $file for writing\n";

# push the standard initial magic, set, and cd to PWD into new output
# array
$pwd = $ENV{"PWD"} 
       or die "make_parexe: Couldn't evaluate variable PWD\n";
push @print_array, "#!/bin/sh\n";
push @print_array, "set -a\n";
push @print_array, "cd $pwd\n";

# Now run though the environment and, in alphabetical order
# set the variables in print_array
foreach $key ( sort keys %ENV ){
  $value = $ENV{$key};

  push @print_array, "$key='$value'\n";
}

# Finally push the LOADMODULE onto the print_array, and write
# out print_array
$loadmodule = $ENV{"LOADMODULE"} 
       or die "make_parexe: Couldn't evaluate variable LOADMODULE\n";
push @print_array, "$loadmodule\n";

print PF @print_array;
close (PF);
