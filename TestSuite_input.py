import sys, getopt
#Basically we have to read a config file and write namelists accordingly
###############################################################################

def CLI_interpreter( argv,
                     config,
                     flags ):

    # get ALL opts & args - where hmfxn are valid options
    # (possible values for "opt") & ":" following valid options
    # represents there is an "arg" to follow
    try:
        #opts, args = getopt.getopt(argv, "hm:f:x:n:", ["ifile=","ofile="])
        opts, args = getopt.getopt(argv, "hm:f:n:")
    # if none of these are present - abort
    except getopt.GetoptError:
        sys.exit(2)

    # Loop over all the options read in, & associate argument
    for opt, arg in opts:

        # if the user wants help (NO arg)
        if opt in ('-h', "--help"):
            #jhan: clean up reportage
            print 'TestSuite.py -f <config_file>\n'
            sys.exit()

        # user config file (arg = <filename>) only valid from CLI
        elif opt in ("-f", "--configfile"):
            config.file.append( arg )
            # reset 1st element indicating CLI arg present
            flags.flag[0] = True
        
    print "CLI requires a config file as an argument."

###############################################################################

def Configfile_interpreter( config,
                            fields,
                            flags,
                            pars ):

    # Declare mutable object => lines read from config file 
    lines = []

    npaths = 0
    # Reads all the lines in the config file specified at CLI
    lines = tuple( open(config.file[0], 'r' ))

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

            # fieldnames
            if lstr.startswith( pars.strchk[0] ):
                cfield = lines[ll].strip().split()
                fields.name.append(cfield[1])

            # paths
            if lstr.startswith( pars.strchk[1] ):
                npaths=npaths+1
                cfield = lines[ll].strip().split()
                fields.path.append(cfield[1])

            # nodes
            if lstr.startswith( pars.strchk[2] ):
                cfield = lines[ll].strip().split()
                config.nodes.append(cfield[1])
                flags.flag[2] = True

            # mapping
            if lstr.startswith( pars.strchk[3] ):
                cfield = lines[ll].strip().split()
                config.map.append(cfield[1])

    # End Loop over read lines[]

    if npaths > 1:
        flags.flag[3] = True
    
    if flags.flag[2] is False:
        print "Aborting. Nodes NOT specified."
        sys.exit()        

    print "Config file Interpreted"

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




