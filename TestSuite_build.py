#!/bin/python

def TestSuite_build:
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
   
   #Start with a clean slate
   cmd = ("/bin/rm -fr " + root + " log" ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/rm -fr /home/599/jxs599/umui_runs/" + UMrun ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/rm -fr /short/p66/jxs599/" + UMsrc ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   cmd = ("/bin/rm -fr /short/p66/jxs599/UM_ROUTDIR/jxs599/" + UMsrc ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
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

   cmd = ("/bin/rm -fr .tmp" ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

   with open(cwd + "/log", "a") as myfile:
      myfile.write("\n\nParallel:\n") 

   # parallel version
   cmd = ("./build_mpi.ksh >> " + cwd + "/log"  ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   cmd = ("/bin/cp cable " + binParallel ) 
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
   cmd = ("/bin/cp -r -p " + cwd + "/" + UMsrc + " /short/p66/jxs599/UM_ROUTDIR/jxs599/ "  ) 
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   os.chdir( "/home/599/jxs599/umui_runs/" + UMrun )
   
   # BUild UM 
   cmd = ("./umuisubmit_compile > " + cwd + "/um_buildlog")
   p = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)
   
   #subprocess.call("ls")

