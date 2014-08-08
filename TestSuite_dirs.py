#!/usr/bin/python
import os

###############################################################################

# Which CABLE applications are you testing (Generally all true)
TestUM = True 
TestMPI = True 
TestSerial = True

###############################################################################

# User specific info (laborious to get from system commnds - implement later)
user = "jxs599"
institute = "599"
project = "p66"

# Users shouldn't need to touch these 
home = "/home/" + institute + "/" + user
short = "/short/" + project + "/" + user

###############################################################################

# For reproducability we use set in stone commited branches

# branch from the CABLE repository branch to test
SVNURL = "https://trac.nci.org.au/svn/cable/trunk"    
# i.e. the root of your CABLE code containing UM,core, offline sub-directories
SVNURLROOT = "/trunk"    

###############################################################################

# Details of the UM job to test with (changes in here will require further playing)
 
# UM Build and Run dirs - this is config dependent on the particular UM run to be tested
# jaaaD is ~ACCESS-1.4 with carbon-cycle switched off 
UMjobID = 'jaaad'
UMjobIDTimeStamp = '-191222944'
UMjobScripts = UMjobID + UMjobIDTimeStamp
UMjobScriptsHome = home + "/umui_runs/" + UMjobScripts
UMrun = short + UMjobID 
UMsrc = short + "/UM_ROUTDIR/" + UMjobID

###############################################################################
# Users shouldn't need to touch anything blow here

# define TestSuite Application dirs

# root of TestSuite Directory
root = os.getcwd()

# root of TestSuiteApps Directory, etc,etc
root_app = root + '/TestSuiteApps'
src = root_app + '/src' 
bin = root_app + '/bin'
run = root_app + '/Run'
cfgs = root + '/TestSuiteConfigs/'
src_root = src + SVNURLROOT
UM = src_root + '/UM'
offline = src_root + '/offline'

#Store for built executables

# Single Sites
binSerial = bin + '/Serial'

# MPI-parallel - generally gswp2
binParallel= bin + '/Parallel'

###############################################################################





