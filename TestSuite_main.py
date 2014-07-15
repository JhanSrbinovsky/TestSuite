#!/usr/bin/python 
__author__ = 'Jhan Srbinovsky'

# Build & Run various CABLE applications described in config file 

# usage: ./TestSuite_main -f <config file>

#import python modules
import sys
import subprocess

#import local, application specific modules
from TestSuite_input import Locate_cfg, Configfile_interpreter

def main(argv):
    
   # insert a break from the CLI for reading
   print "\n"
   
###############################################################################
   # Locate config file (ifile), log file (ofile)
   ifile = [] 
   ofile = []
   Locate_cfg( argv, ifile, ofile )
   #print "ifile ", ifile[0] 
   #print "ofile ", ofile[0] 
   
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
   #for i in range( 1 ):   
   #   cfg.append( Configs() )
   cfg = Configs() 
   
   # Read config file and fill cfg. class
   Configfile_interpreter( ifile, cfg )
   
   # Print what we are going to do cfg. class
   print "Applications to be evaluated:" 
   for i in range( len( cfg.name ) ):   
      print "\n",cfg.name[i] 

   # print class all at once 
   #attrs = vars( cfg )
   # now dump this in some way or another
   #print ', '.join("%s: %s" % item for item in attrs.items())
   sys.exit()
 
###############################################################################
    
   
   #execute system commands
   print "In current format, dumped files are longitudinally complete per node."
   print "Therefore cannot be processed by same methodi used here."
   cmd1 = ("/bin/cp " + pars.path[0] + "/longitude001.bin " 
          + pars.path[1] + "/longitude.bin")
   cmd2 = ("/bin/cp "+ pars.path[0] + "/longitude001.dat " 
          + pars.path[1] + "/longitude.dat")
   
   #subprocess.call(cmd)
   p = subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True)
   (output, err) = p.communicate()
   p = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)
   (output, err) = p.communicate()
   
   # Get user input from both CLI and config file
   user_input( argv,
               config,
               pars,
               fields,
               flags )
   
   # Operate on data
   driver(     fields,
               pars,
               config,
               flags )
   
   
   
   
   
   





# comment out if interactive
#################################################################################
if __name__ == "__main__":
   main(sys.argv[1:])

################################################################################


