#!/usr/bin/python

#import python modules
import os
import sys
import subprocess

#import local, application specific modules
from TestSuite_dirs import root, root_app, src, bin, run, binSerial, binParallel 
from TestSuite_dirs import trunk, UM, offline, UMrun, UMsrc, SVNURL, cfgs
from TestSuite_dirs import TestUM, TestSerial, TestMPI 

def TestSuite_runner( cfg ):

###############################################################################

   os.chdir( root )
   # Run Applications

   for i in range( len( cfg.path ) ):  
      #For serial runs      
      if(str(cfg.mode[i]) == '1'):

         if TestSerial is True: 
            # cp executable and namelist to rundir
            rundir = str( run + "/" + cfg.path[i] )

            cmd = ("/bin/cp " + binSerial + "/cable " + rundir ) 
            p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
            cmd = ("/bin/cp " + cfgs + cfg.path[i] + "/cable.nml " + rundir ) 
            p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

            # GoTo rundir and execute
            os.chdir( rundir )
            cmd = ("./cable" ) 
            p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
         
      #For parallel runs      
      if(str(cfg.mode[i]) == '2'):
         
         if TestMPI is True: 
            # cp executable and namelist to rundir
            rundir = str( run + "/" + cfg.path[i] )

            cmd = ("/bin/cp " + binParallel + "/cable-mpi " + rundir ) 
            p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
            cmd = ("/bin/cp run_cable-mpi " + rundir  ) 
            p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
            #cmd = ("/bin/sed \'s/RUNDIR/" + rundir + "/ \'" + cfgs + cfg.path[i] + "/run_cable-mpi" ) 
            #p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
            cmd = ("/bin/cp " + cfgs + cfg.path[i] + "/cable.nml " + rundir ) 
            p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
        
            # GoTo rundir and execute
            os.chdir( rundir )
            cmd = ("qsub run_cable-mpi" ) 
            p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
 
      #For UM runs      
      # UM executable is already in the right place         
      if(str(cfg.mode[i]) == '3'):
         os.chdir( "/home/599/jxs599/umui_runs/" + UMrun )
         
         if TestUM is True: 
            # Run UM 
            cmd = ("qsub umuisubmit_run > " + root + "/um_runlog")
            p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
               



