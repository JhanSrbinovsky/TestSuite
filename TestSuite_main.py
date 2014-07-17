#!/usr/bin/python 
__author__ = 'Jhan Srbinovsky'

# Build & Run various CABLE applications described in config file 

# usage: ./TestSuite_main -f <config file>

#import python modules
import os
import sys
import subprocess

#import local, application specific modules
from TestSuite_input import Locate_cfg, Configfile_interpreter
from TestSuite_build import TestSuite_builder, CleanSlate 
from TestSuite_run import TestSuite_runner 
#from TestSuite_run import TestSuite_run 

def main(argv):

   # insert a break from the CLI for reading
   print "\n"
   
###############################################################################
   
   # Locate config file (ifile), log file (ofile)
   # Defaults to using default.cfg as input, default.out as output (currently no output) 
   ifile = [] 
   ofile = []
   Locate_cfg( argv, ifile, ofile )
   
###############################################################################

   # Read config file
   
   # class to hold all config file input
   class Configs(object):
      def __init__(self):
         # Declare flags as mutable:
         self.mode = []
         self.name = []
         self.path = []
   cfg = []
   cfg = Configs() 
   
   # Read config file and fill cfg. class
   Configfile_interpreter( ifile, cfg )
   
   # Print what we are going to do cfg. class
   print "Applications to be evaluated:" 
   for i in range( len( cfg.name ) ):   
      print "\n",cfg.name[i] 
 
###############################################################################
   
   # Execute Applications 
   
   print "\nSet up structure for building/running all Apps described in " + \
   "config file. These will be cleaned up following execution of the TestSuite"

   #Start with a clean slate
   CleanSlate()

   # Build Applications
   TestSuite_builder( cfg )
   
   # Run Applications
   TestSuite_runner( cfg )

# we may simply be able to use run.ksh here
# worry about qsub later. first run on a single node   
  
## here it may be easier to use ksh scripts   
#   for i in range( len( cfg.name ) ):   
#      cmd = ( "./TestSuite.ksh " +  cfg.path[i] ) 
#      p = subprocess.check_callPopen(cmd, stdout=subprocess.PIPE, shell=True)
#     
# 


###############################################################################

   # Evaluate Success of Applications and write to ofile [default.out] 
   
   
   
   
#Atm_Step: Timestep                      2
   
   





# comment out if interactive
#################################################################################
if __name__ == "__main__":
   main(sys.argv[1:])

################################################################################


