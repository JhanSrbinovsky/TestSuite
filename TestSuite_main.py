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
   bin = root + '/bin'
   run = root + '/Run'
   binSerial = bin + '/Serial'
   binParallel= bin + '/Parallel'
   binUM = bin + '/UM'
   trunk = src + '/trunk'
   UM = trunk + '/UM'
   offline = trunk + '/offline'
   UMrun = 'jaaad-191222944'
   UMsrc = 'jaaad'
    
###############################################################################
   # mkdir structure
   print "mkdir structure\n"
   cmd = ("/bin/mkdir -p " + root ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + src ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + bin ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + binSerial ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + binParallel ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + binUM ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + run ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   for i in range( len( cfg.path ) ):   
      cmd = ("/bin/mkdir -p " + run + "/" + cfg.path[i] ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
               
###############################################################################

   # checkout trunk (or URL to test - remove hardwiring)
   os.chdir(src)
   #use this on raijin
   cmd = ("/usr/bin/svn co " + SVNURL ) 
   #jiggle
   #cmd = ("/opt/subversion/bin/svn co " + SVNURL ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
###############################################################################

   with open(cwd + "/log", "a") as myfile:
      myfile.write("BUILDS:\n\n") 
      myfile.write("Serial:\n") 

###############################################################################

   # build models
   print "build model\n"

   # overwrite build scripts checked out 
   cmd = ("/bin/cp " + cwd + "/build.ksh " + offline ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/cp " + cwd + "/build_mpi.ksh " + offline ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   # offline
   os.chdir(offline)
   
   # serial version
   #load netcdf in build script 
   cmd = ("./build.ksh >> " + cwd + "/log" ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

   cmd = ("/bin/cp cable " + binSerial ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

   with open(cwd + "/log", "a") as myfile:
      myfile.write("\n\nParallel:\n") 

   # parallel version
   cmd = ("./build_mpi.ksh >> " + cwd + "/log"  ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   # UM 
   with open(cwd + "/log", "a") as myfile:
      myfile.write("Build libcable first....\n\n") 
   os.chdir(UM)
   cmd = ("./build.ksh >> " + cwd + "/log" ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   with open(cwd + "/log", "a") as myfile:
      myfile.write("Then move to UM build....\n\n") 
   
   # cp UM runscripts to execute from
   cmd = ("/bin/cp -r " + cwd + "/" + UMrun + " /home/599/jxs599/umui_runs/ " ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   # cp UM Extracted src directory 
   cmd = ("/bin/cp -r -p" + cwd + "/" + UMsrc + " /short/p66/jxs599/UM_ROUTDIR/jxs599/ "  ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   os.chdir( "/home/599/jxs599/umui_runs/" + UMrun )
   
   # BUild UM 
   cmd = ("./umuisubmit_compile > um_buildlog")
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   
   
   subprocess.call("ls")

   sys.exit()
       
###############################################################################
   
   # Run Applications

   for i in range( len( cfg.path ) ):   
      if(str(cfg.mode[i]) == '1'):
         cmd = ("/bin/cp " + binSerial + "/cable " + run + "/" + cfg.path[i] ) 
         p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
               
   # once we have the lib we have to build the UM as well 


   #run models
#we may simply be able to use run.ksh here
# worry about qsub later. first run on a single node   
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


