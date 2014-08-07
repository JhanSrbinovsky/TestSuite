#!/usr/bin/python

#import python modules
import os
import sys
import subprocess

#import local, application specific modules
from TestSuite_dirs import root, root_app, src, bin, run, binSerial, binParallel 
from TestSuite_dirs import trunk, UM, offline, UMrun, UMsrc, SVNURL, cfgs
from TestSuite_dirs import TestUM, TestSerial, TestMPI 

def TestSuite_builder( cfg ):

###############################################################################

   # mkdir structure
   print "mkdir structure\n"
   cmd = ("/bin/mkdir -p " + root_app ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + src ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + bin ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + binSerial ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + binParallel ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/mkdir -p " + run ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   for i in range( len( cfg.path ) ):   
      cmd = ("/bin/mkdir -p " + run + "/" + cfg.path[i] ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
               
################################################################################

   # checkout trunk (or URL to test - remove hardwiring)
   os.chdir(src)
   #use this on raijin
   cmd = ("/usr/bin/svn co " + SVNURL ) 
   #jiggle
   #cmd = ("/opt/subversion/bin/svn co " + SVNURL ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
###############################################################################

   with open(root + "/log", "a") as myfile:
      myfile.write("BUILDS:\n\n") 
      myfile.write("Serial:\n") 

###############################################################################

   # build models

   # overwrite build scripts checked out 
   cmd = ("/bin/cp " + root + "/build.ksh " + offline ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/cp " + root + "/build_mpi.ksh " + offline ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   # serial version
   if TestSerial is True: 
      # offline
      os.chdir(offline)
   
      print "build serial model\n"
      cmd = ("./build.ksh >> " + root + "/log" ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

      print "build serial model finished\n"
      
      cmd = ("/bin/cp cable " + binSerial ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

      print "copied serial model to the bin/ directory\n"
      
      cmd = ("/bin/rm -fr .tmp" ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

      with open(root + "/log", "a") as myfile:
         myfile.write("\n\nParallel:\n") 

   # parallel version
   if TestMPI is True:
      # offline
      os.chdir(offline)
   
      print "build model\n"
      cmd = ("./build_mpi.ksh >> " + root + "/log"  ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
      
      print "build parallel model finished\n"
      
      cmd = ("/bin/cp cable-mpi " + binParallel ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

      print "copied parallel model to the bin/ directory\n"

      # UM 
   if TestUM is True: 
      with open(root + "/log", "a") as myfile:
         myfile.write("Build libcable first....\n\n") 
      os.chdir(UM)
      cmd = ("./build.ksh >> " + root + "/log" ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
      
      with open(root + "/log", "a") as myfile:
         myfile.write("Then move to UM build....\n\n") 
   
      # cp UM runscripts to execute from
      cmd = ("/bin/cp -r " + root + "/" + UMrun + " /home/599/jxs599/umui_runs/ " ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
      
      # cp UM Extracted src directory 
      cmd = ("/bin/cp -r -p " + root + "/" + UMsrc + " /short/p66/jxs599/UM_ROUTDIR/jxs599/ "  ) 
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
      
      os.chdir( "/home/599/jxs599/umui_runs/" + UMrun )
      
      # BUild UM 
      cmd = ("./umuisubmit_compile > " + root + "/um_buildlog")
      p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
################################################################################

def CleanSlate():   
   cmd = ("/bin/rm -fr " + root_app + " log" ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/rm -fr /home/599/jxs599/umui_runs/" + UMrun ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/rm -fr /short/p66/jxs599/" + UMsrc ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/rm -fr /short/p66/jxs599/UM_ROUTDIR/jxs599/" + UMsrc ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/rm -f " + root + "/*log*" ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

################################################################################



