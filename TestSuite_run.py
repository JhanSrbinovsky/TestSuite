#!/usr/bin/python

#import python modules
import os
import sys
import subprocess

#import local, application specific modules
from TestSuite_dirs import cwd, root, src, bin, run, binSerial, binParallel 
from TestSuite_dirs import binUM, trunk, UM, offline, UMrun, UMsrc, SVNURL

def TestSuite_runner( cfg ):

###############################################################################

   os.chdir( cwd )
   # Run Applications

   for i in range( len( cfg.path ) ):  
      #For serial runs      
      if(str(cfg.mode[i]) == '1'):

         # cp executable and namelist to rundir
         rundir = str( run + "/" + cfg.path[i] )

         cmd = ("/bin/cp " + binSerial + "/cable " + rundir ) 
         p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
         cmd = ("/bin/cp " + "TestSuiteConfigs/" + cfg.path[i] + "/cable.nml " + rundir ) 
         p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

         # GoTo rundir and execute
         os.chdir( rundir )
         cmd = ("./cable" ) 
         p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
         
      #For parallel runs      
      if(str(cfg.mode[i]) == '2'):
         cmd = ("/bin/cp " + binParallel + "/cable " + run + "/" + cfg.path[i] ) 
         p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

      #For UM runs      
      # UM executable is already in the right place         
      if(str(cfg.mode[i]) == '3'):
         os.chdir( "/home/599/jxs599/umui_runs/" + UMrun )
         
         # Run UM 
         cmd = ("qsub umuisubmit_run > " + cwd + "/um_runlog")
         p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
               



