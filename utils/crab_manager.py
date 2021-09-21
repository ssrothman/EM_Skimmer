# -------------------------------------
# Script to generate CRAB jobs
# -------------------------------------

#!/usr/bin/env python

import os, sys

from yaml_utils import *
from collections import OrderedDict

from datetime import datetime
from time_utils import *


class crab_job_manager():
  
    def __init__(self, file_ = ''):

        self.config_file = file_
        self.config = None
        self.valid_config = False
        self.parameter_list = ['InputDataset', 'OutputFile','RequestName','WorkArea',\
                 'TransferOutputs','TransferLogs','PsetName','MaxJobRunTime',\
                 'InputDBS','Splitting','UnitsPerJob','FilePath','StorageSite',\
                 'UseWhitelist', 'UseBlacklist',\
                 'Whitelist','Blacklist', 'IgnoreGlobalBlacklist']

    def check_config(self):

        for key_ in self.parameter_list:
            if key_ not in self.config:
                print 'Specify '+key_+' in the config file!'
                return False

        return True

    def read_config_file(self, file_):

        self.config_file = file_
        self.config = ordered_load(open(file_), Loader=yaml.SafeLoader)
        if not self.check_config():
            print 'Loading ocnfig failed. Check example_config.yml for help.'
            self.valid_config = False

        else: self.valid_config = True

    def produce_crab_submission_script(self):

        if not self.valid_config:
            print 'Please load config file first!'
            return

        script_name = 'crab_{}.py'.format(self.config['RequestName'])

        with open('crab_template.txt', 'read') as tmp:
            with open(script_name, 'w') as sub:

                plain_text = tmp.readlines()
                new_text = []

                for line in plain_text:
                    newline = line
                    if 'Template file for CRAB job submission' in line:
                        newline = '# CRAB Job: {}\n'.format(print_creation_timestamp())
                    for key_ in self.parameter_list:
                        par_ = self.config[key_]
                        if '<{}>'.format(key_) in line:
                            # Check for the type of the parameter
                            # If the parameter is a string replace it with a string
                            # Else if it's a float or an integer or a list, convert to a string
                            if type(par_)==str: newline = line.replace('<{}>'.format(key_), "'{}'".format(par_))
                            else:
                                newline = line.replace('<{}>'.format(key_), str(par_))
                                if key_=='Blacklist' and self.config['UseBlacklist']==False: newline = '#'+newline
                                if key_=='Whitelist' and self.config['UseWhitelist']==False: newline = '#'+newline
                    new_text.append(newline)


                for line in new_text:
                    sub.write(line)

        return script_name


    def setup_crab(self, file_):
  
        self.read_config_file(file_)
        _ = self.produce_crab_submission_script()

        return _
