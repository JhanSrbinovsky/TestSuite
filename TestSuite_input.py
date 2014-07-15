#/usr/bin/python 
import sys, getopt

#Basically we have to read a config file and write namelists accordingly
###############################################################################

def Locate_cfg( argv, ifile, ofile ):
   inputfile = 'default.cfg'
   outputfile = 'default.log'
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'test.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   
   ifile.append( inputfile )
   ofile.append( outputfile )

###############################################################################

def Configfile_interpreter( ifile, cfg ):

   # Declare mutable object => lines read from config file 
   lines = []
   
   # Reads all the lines in the config file specified at CLI
   lines = tuple( open(ifile[0], 'r' ))
   
   # Loop over all of these read lines
   for ll in range( len(lines) ):
   
      # Create a string from each line
      lstr = str(lines[ll])
      # nullify case sensitivity in that string
      lstr = lstr.lower()
      
      # if it is meant to be interpreted - interpret based on what
      # it starts with. See configfile example format
      # check if the string is a comment
      if not lstr.startswith('#'):
      
         # Mode
         if lstr.startswith( "mode" ):
            cfield = lines[ll].strip().split()
            cfg.mode.append(cfield[1])
         
         # Name
         if lstr.startswith( "name" ):
            cfield = lines[ll].strip().split()
            cstr = ""
            for istr in range( 1,len(cfield) ):
               cstr = cstr + cfield[istr] 
            cfg.name.append(cstr)
         
         # Path
         if lstr.startswith( "path" ):
            cfield = lines[ll].strip().split()
            cfg.mode.append(cfield[1])
      
      # End IF not comment line

   # End Loop over read lines[]

###############################################################################

