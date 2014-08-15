#!/usr/bin/python

from optparse import OptionParser



#############################################
#
#   command line argument parser
#

class CLParser:


    def __init__(self):
        add_help_option = "./assist -a ACTION [optional parameters]"
        
        self.parser = OptionParser(add_help_option, conflict_handler='resolve')
        
        self.parser.add_option("-a", "--action-name", dest="action", default=None,
                          help="Action name", metavar="ACTIONNAME")
        
        self.parser.add_option("-D", "--debug-level", dest="debug", default=None,
                          help="Debug level", metavar="DEBUGLEVEL")
        
        

    def parse_args(self):
        (self.options, self.args) = self.parser.parse_args()



    def get_options(self):
        
        return self.options


    def print_usage(self):
        self.parser.print_usage()

        
    
#
##### end of command line parser
