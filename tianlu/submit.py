""" This module contains the base submit class for submitting jobs to
the batch.  Classes that inherit from this will need to overwrite the
self.loop() method to write the self.prog_dict for job submission.
self.prog_dict[key] = value has key = job_name, value = command
"""

import os
import time
import subprocess
import shlex
import getpass
from argparse import ArgumentParser
from tianlu import utils


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


class submit(object):

    def __init__(self):
        # Variables that can be overwritten
        self.jobdir = 'jobs/'
        # allows us to overwrite job output filepath
        # default is jobdir/job_str_time.log
        self.outputarg = None
        # self.prog_dict is a dictionary of the form:
        # {'job_arg':'program_to_run'}
        self.prog_dict = {}

        self._cmd = 'Qsub -e -l lnxfarm -c TMPDIR'
        self._homedir = os.environ['HOME']

        # parse options
        # we don't call self.parser.parse_args() in __init__ in case the
        # child class of submit needs to define additional options.  This
        # gives more flexibility in writing child classes and allows the
        # user to call parse_args() when all options have been defined.
        self.parser = ArgumentParser()
        self.parser.add_argument('-n',
                               dest='dry_run',
                               action='store_true',
                               default=False,
                               help='dry run')
        self.parser.add_argument('-m',
                               dest='email',
                               type=str,
                               default=None,
                               help='send email when job ends')


    def add_job(self, cmd, job_name=None):
        """ Add an individual job to the job_dict.  This way it's not
        necessary to subclass this class and overwrite loop().
        """
        if job_name is None:
            job_name = os.path.splitext(os.path.basename(cmd))[0]

        self.prog_dict[job_name] = cmd


    def loop(self):
        """ Method to be overwritten.  Generates self.prog_dict
        """
        pass


    def run(self):
        (options, args) = self.parser.parse_known_args()
        if options.email is None:
            mailarg = ''
        else:
            # see qsub man page for details.  a(abort), b(beginning), e(end)
            mailarg = ' -m ae -M ' + options.email

        self.loop()
        print 'NJOBS:', len(self.prog_dict)
        for job_str, aProg in self.prog_dict.iteritems():
            if not self.outputarg:
                output = ' -o ' + self.jobdir + job_str + \
                    '_' + str(int(time.time())) + '.log '
            else:
                output = self.outputarg
                self.jobdir = os.path.dirname(shlex.split(output)[1])
            jobarg = ' -N ' + job_str
            oscmd = self._cmd + mailarg + jobarg + output + aProg
            print oscmd
            if not options.dry_run:
                utils.make_dirs_if_needed(self.jobdir)
                utils.queue_check(getpass.getuser())
                subprocess.call(shlex.split(oscmd))
                time.sleep(5)
