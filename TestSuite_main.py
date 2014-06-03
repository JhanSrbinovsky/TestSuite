#!/usr/bin/python
__author__ = 'Jhan Srbinovsky'

# Process binary fields dumped per node from CABLE into global netcdf files.

# usage: ./bins2ncdf_main -f <config file>

# Required format for config file described in example configfile

#import python modules
import sys
import subprocess

#import local, application specific modules
from bins2ncdf_input import user_input 
from bins2ncdf_driver import driver 

def main(argv):
    
    # insert a break from the CLI for reading
    print "\n"
    
    # Defaults defined. Ability to overwrite these is depracated. Therefore 
    # change here if necessary. Default directories are given relative to ./ 

    class Params(object):
        def __init__(self):

            # Check per line "label" in configfile defining the type
            # of the following argument 
            self.strchk =   [ "FieldName",
                              "Path",
                              "Nodes",
                              "Mapping",
                              "MapPath" ]

            # Directories: mapping data, cocatenated data, exececutables 
            self.path =     [ "data/mapping",
                              "data/catted",
                              "bin/cat_Nnodes", 
                              "bin/ncdf_main", 
                              "data/ncdf" ] 

            # Fieldnames of mapping data 
            self.mapfield = [ "latitude",
                              "lat_index",
                              "lon_index",
                              "tile_frac",
                              "tile_index" ]

    # Declare flags as mutable:
    # [0] - configfile, [1] - mapping, [2] - nodes, [3] - N paths 
    class Flags(object):
        def __init__(self):
            self.flag = []

    class Fields(object):
        def __init__(self):
            self.name  = []
            self.path  = []

    
    # Can be set either via CLI or commandline: trumped by configfile
    class Configs(object):
        def __init__(self):
            self.file = []      # Can be set only via CLI 
            self.map  = []
            self.nodes  = []

    pars = Params()
    flags = Flags()
    fields = Fields()
    config = Configs() 

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











################################################################################
if __name__ == "__main__":
   main(sys.argv[1:])

################################################################################


