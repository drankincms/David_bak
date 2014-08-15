#!/usr/bin/env python


#####################################################
#
#   batch job and other actions manager
#
#   routines for submitting and managing CRAB,
#   Condor jobs and related activities
#


#
# singleton base class for all actions
# serves as the action manager
#
class assistAction:
    theInstance = None

    @classmethod
    def instance(cls):
        if cls == assistAction:
            if assistAction.theInstance == None:
                assistAction.theInstance = assistAction()

        return assistAction.theInstance
    

    def __init__(self):

        if self.__class__ == assistAction:
            if assistAction.theInstance:
                raise singletonException('The assistAction class is a singleton. Use assistAction.instance() to get its instance')

            assistAction.theInstance = self

        self.name = 'assistAction'
        self.legend = '[Action manager]:'

        print self.legend, 'loaded'


        if self.__class__ == assistAction:
            self.action_dict = {}
        
    def register(self):
        print self.legend, 'registering', self.__class__.__name__
        assistAction.instance().action_dict[self.get_name()]=self

    def action(self, ws = None):
        print self.legend, 'action'

    def get_legend(self):
        return self.legend

    def get_name(self):
        return self.name

    def print_actions(self):
        print self.legend, 'Available actions:'
        for act in self.action_dict:
            print self.legend, '  ', self.action_dict[act].get_name()

    def action_names(self):
        return self.action_dict.keys()
