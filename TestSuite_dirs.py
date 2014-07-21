#!/usr/bin/python
import os

###############################################################################
TestUM = False
TestMPI = False
TestSerial = True 
###############################################################################

# define dirs

# branch from the CABLE repository branch to test
SVNURL = "https://trac.nci.org.au/svn/cable/trunk"    

# root of TestSuite Directory
root = os.getcwd()

# root of TestSuiteApps Directory, etc,etc
root_app = root + '/TestSuiteApps'
src = root + '/src' 
bin = root + '/bin'
run = root + '/Run'
trunk = src + '/trunk'
UM = trunk + '/UM'
offline = trunk + '/offline'

#Store for built executables

# Single Sites
binSerial = bin + '/Serial'

# MPI-parallel - generally gswp2
binParallel= bin + '/Parallel'

# UM Build and Run dirs - this is config dependent on the particular UM run to be tested
# jaaaD is ~ACCESS-1.4 with carbon-cycle switched off 
UMrun = 'jaaad-191222944'
UMsrc = 'jaaad'

###############################################################################




