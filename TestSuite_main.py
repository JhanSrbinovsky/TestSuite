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

def main(argv):
   SVNURL = "https://trac.nci.org.au/svn/cable/trunk"    
   # insert a break from the CLI for reading
   print "\n"
   
###############################################################################
   # Locate config file (ifile), log file (ofile)
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

   cwd = os.getcwd()
   root = cwd + '/TestSuiteApps'
   src = root + '/src' 
   trunk = src + '/trunk'
   UM = trunk + '/UM'
   offline = trunk + '/offline'

   # mkdir structure
   print "mkdir structure\n"
   cmd = ("/bin/mkdir -p " + root ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   #p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + src ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   #p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
   
   # checkout trunk (or URL to test - remove hardwiring)
   os.chdir(src)
   #use this on raijin
   #cmd = ("/usr/bin/svn co " + SVNURL ) 
   #jiggle
   cmd = ("/opt/subversion/bin/svn co " + SVNURL ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   #p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
   
   # build models
   print "build model\n"
   cmd = ("/bin/cp " + cwd + "/build.ksh " + offline ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   # offline
   os.chdir(offline)
   
   # serial version
   print "\n" 
   print "\n" 
   cmd = ("./build.ksh > " + cwd + "/log" ) 
   print "\n" 
   print "\n" 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   subprocess.call("ls")
   
   print "\n" 
   print os.getcwd()
   sys.exit()     
   
   #p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
   
   ## parallel version
   #cmd = ("./build_mpi.ksh" ) 
   #p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   ##p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
   #
   ## UM 
   #os.chdir(UM)
   #cmd = ("./build.ksh" ) 
   #p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   ##p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
  
   # once we have the lib we have to build the UM as well 

   sys.exit()     
  
  
#   for i in range( len( cfg.name ) ):   
#      cmd = ( "./TestSuite.ksh " +  cfg.path[i] ) 
#      p = subprocess.check_callPopen(cmd, stdout=subprocess.PIPE, shell=True)
#     
# 
## here it may be easier to use ksh scripts   

   
   
   
   
   





# comment out if interactive
#################################################################################
if __name__ == "__main__":
   main(sys.argv[1:])

################################################################################


