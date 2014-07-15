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
   #print 'Input file is "', inputfile
   #print 'Output file is "', outputfile

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

def user_input( argv,
                config,
                pars,
                fields,
                flags ):

    # Nullify case sensitivity
    pars.strchk[0] = pars.strchk[0].lower()
    pars.strchk[1] = pars.strchk[1].lower()
    pars.strchk[2] = pars.strchk[2].lower()
    pars.strchk[3] = pars.strchk[3].lower()
    pars.strchk[4] = pars.strchk[4].lower()

    # initialize first element of object
    for i in range(4):
       flags.flag.append( False )

    # CALL function to read from command line
    CLI_interpreter( argv,
                     config,
                     flags )

    if flags.flag[0] is False:
        print "\nA config file must be specified.\n "
        print "\t .MapMe.py -f <filename>\n"
        print "There is no other mechanism to describe what fields you want processed."
        if flags.flag[1] is False:
            print "\nShould mapping be performed? Unless specified at the CLI, \n"
            print "\t .MapMe.py -m <T/F>\n"
            print "it has to be in your config file. Otherwise it is assumed True.\n "
        if flags.flag[2] is False:
            print "Unless specified at the CLI, \n"
            print "\t .MapMe.py -n <nodes>\n"
            print "this has to be in your config file. Otherwise we abort.\n "

    else:
        print 'Opening config file - ', config.file[0]
        # func: read config file
        Configfile_interpreter( config,
                                fields,
                                flags,
                                pars )

        if flags.flag[2] is False:
            print "We have no idea how many nodes you used. " + \
                  "Tell us via the CLI or in your config file?\n"
            sys.exit()

        print "\nOK... to summarize what we are about to work on, " + \
              "as gathered from your CLI args and config file: \n"
        print 'Originally data was dumped from ', str(config.nodes[0]), ' nodes.'

        if config.map[len(config.map)-1] is True:
            print 'We are going to process the mapping data in path:'
            print '\t',pars.path[0]
        else:
            print 'We are NOT going to process the mapping data'

        print 'We will map the following fields: '
        for ff in range( len(fields.name)):
            print '\t', fields.name[ff]
        if flags.flag[3] is True:
            print '\n Corresponding to paths'
            for pp in range( len(fields.path)):
                print '\t', fields.path[pp]
        else:
            print "\nAll fields are in: "
            print '\t', fields.path[0]




