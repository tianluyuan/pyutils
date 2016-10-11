""" This is a subclass of submit that handles job submission for processing
nd280 data.
"""
import shutil
import time
import os
import shlex
import inspect
from natsort import natsorted
from tianlu import utils, directories
from tianlu.submit import submit
from tianlu.analysis import highland_configurations as hc
from tianlu.analysis import utils as tautils


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


# Group executables and scripts together by arguments, i.e. all
# scripts with the same number and format of args should be in a
# single list below.  Each list should have its own _add_ function in
# CmdBuilder
#
# All highland/highland2 analyses should go here
HIGHLAND_SCRIPTS = ['RunNumuCCAnalysis',
                    'RunP0DNumuCCAnalysis',
                    'RunP0DNumuCCQEAnalysis']

# All scripts dependent on runAnalyses.py
RUNANALYSES_SCRIPTS = ['GetEvents',
                       'GetCosmicsEvents',
                       'GetSandEvents',
                       'GetCosmicsInBeam']

# All oaEvent scripts written to process reco files
TIANND280_SCRIPTS = ['getMCMTimes']

# Individual script in xsTool for creating T2KReWeight weights
CREATE_WEIGHTS_SCRIPTS = ['create_weights_file']


class CmdBuilder(object):
    """Class for building the cmdline program for individual analysis
    scripts.  We need this since different scripts will have different
    cmdline args, input/output filepaths, etc.  Each unique cmdline
    format must be created by a method that begins with '_add_' and
    updates the self._prog_dict variable to contain the script name as
    key and the job template as the dict's value
    """
    def __init__(self):
        """ Create CmdArgs object for every class of scripts that have different
        commandline arguments.  This will include the prog_cmd with various args
        passed to the program, the input-specifier which contains the path to
        the sorted subdir, and the output filepath.

        Each subdir has an associated .txt file so we can
        specify which input file to use with the final string (e.g. '/*.root' for
        runAnalyses.py or '.txt' for using input filelists with highland).
        """
        # This is unique to each job and common to the input, output, and
        # filenaming
        self._unique_fmt = '{prod}/{dtype}/{gen}/{run}'
        self._outfile_fmt = self._unique_fmt.replace('/', '_')

        # path templates
        self.data_dir_template = os.path.join('{base}',
                                              self._unique_fmt,
                                              '{ftype}')

        outfile_template = '{script}_'+ self._outfile_fmt+'_{subdir}.root'
        self._output_common = os.path.join('{script}',
                                           self._unique_fmt,
                                           outfile_template)
        self._output_weights = os.path.join('{microtree}',
                                            self._unique_fmt,
                                            '{subdir}/')

        # analysis program template (not used by create_weights)
        self._pretemp = ('{prog_cmd} {input_specify} -o {outfile} && '
                         'CUStageOut {outfile} {output}')

        # create and setup _prog_dict
        self._prog_dict = {}

        # loop over '_add_' members and execute all of them
        # uses reflection
        for name, member in inspect.getmembers(self,
                                               predicate=inspect.ismethod):
            if name.startswith('_add_'):
                member()


    def _update_prog_dict(self, template, scripts):
        for s in scripts:
            self._prog_dict[s] = template


    def _add_highland_templates(self):
        """ Create command templates for each program in HIGHLAND_SCRIPTS and
        fill into self._prog_dict.

        Since there is now highland and highland2, and since their programs
        have the same name but we want to save outputs into different
        directories, we must check the highland version for each individual
        program.
        """
        for script in HIGHLAND_SCRIPTS:
            highland_version = tautils.exe_highland_version(script+'.exe')
            if highland_version == 2:
                output = os.path.join(directories.HIGHLAND2_STORE,
                                      self._output_common)
            elif highland_version == 1:
                output = os.path.join(directories.HIGHLAND_STORE,
                                      self._output_common)

            template = self._pretemp.format(prog_cmd='{script}.exe',
                                            input_specify='{inp}.txt -p {paramfile}',
                                            outfile=os.path.basename(output),
                                            output=output)
            self._update_prog_dict(template, [script])


    def _add_runanalyses_templates(self):
        output = os.path.join(directories.RUNANALYSES_STORE,
                              self._output_common)
        prog = 'cp $SGE_O_WORKDIR/{script}.py ./ && runAnalyses.py --batch -a {script}'

        template = self._pretemp.format(prog_cmd=prog,
                                        input_specify='{inp}/*.root',
                                        outfile=os.path.basename(output),
                                        output=output)
        self._update_prog_dict(template, RUNANALYSES_SCRIPTS)


    def _add_tiannd280_templates(self):
        output = os.path.join(directories.TIANND280_STORE,
                              self._output_common)

        template = self._pretemp.format(prog_cmd='getMCMTimes.exe -g -t C',
                                        input_specify='{inp}/*.root',
                                        outfile=os.path.basename(output),
                                        output=output)
        self._update_prog_dict(template, TIANND280_SCRIPTS)


    def _add_create_weights_templates(self):
        """ This creates cmd_temps that will run create_weights_file.sh in the
        xsTool.  It requires a few more arguments than the other scripts.

        Copies script over to TMPDIR and outputs everything locally, then
        uses CUStageOut to move files into permanent storage area (outdir).
        """
        pretemp = ('create_weights_file_wrapper.sh {work_dir} {inputtemp} '
                   '{banfffile} {oaanallist} {massfile} {outdir}')

        inp_store = self._get_highland_store()
        inputtemp = os.path.join(inp_store,
                                 self._output_common.replace('{script}',
                                                             '{microtree}'))

        banfffile = hc.BANFF_FILE
        massfile = hc.MASS_FILE
        outdir = os.path.join('{weights_store}', self._output_weights)
        oaanallist = '{inp}.txt'
        if os.environ.has_key('XSTOOLROOT'):
            work_dir = os.path.join(os.getenv('XSTOOLROOT'), 'work')
        else:
            print '$XSTOOLROOT not set!'
            print 'Not setting up cmd dict for create_weights_file!'
            return

        template = pretemp.format(work_dir=work_dir,
                                  inputtemp=inputtemp,
                                  banfffile=banfffile,
                                  massfile=massfile,
                                  outdir=outdir,
                                  oaanallist=oaanallist)
        self._update_prog_dict(template, CREATE_WEIGHTS_SCRIPTS)


    def _get_highland_store(self):
        """Gets the input and output storage areas based on whether the
        environment is setup for processing highland2 files or
        highland1

        Returns inp_store (path to highland store), out_store (path to
        weights store)
        """
        if tautils.env_is_highland2():
            inp_store = directories.HIGHLAND2_STORE
        else:
            inp_store = directories.HIGHLAND_STORE

        return inp_store


    def build(self, job):
        full_template = self._prog_dict[job['script']]

        cmd = full_template.format(**job)
        return cmd


class SubmitAnalysis(submit):
    """ This class basically concatenates a bunch of
    strings which are then passed to the os to be
    submitted as commands to the cluster.

    Description of variables and their requirements
    1) script_list: List of python scripts in cwd
    that will be used to process oaAnalyses files
    and skim events based on some base cut.  Should
    typically create a skimmed outFile with relevant,
    useful variables
    2) data_list: Can be 'rdp' and/or 'mcp'.  Tells
    submit.py to process rdp and/or mcp files.
    3) generator_list: Can be 'beam' (for rdp beam)
    'neut'/'genie' (for mcp beam) or 'cosmics'
    (for both rdp and mcp cosmics).
    4) run_list: For neutrino events, must be of the
    form 'RUN1', 'RUN2water' etc.  For cosmic
    files, must be of the form '201002water' (specifying
    the nd280 geometry used for mc and recon).
    """
    def __init__(self, config_dict):
        super(SubmitAnalysis, self).__init__()
        self.parser.add_argument('--dev',
                               dest='dev',
                               action='store_true',
                               default=False,
                               help='process a small test sample')
        self.parser.add_argument('-f',
                               dest='force',
                               action='store_true',
                               default=False,
                               help='suppress prompt and force processing')

        self._use_hadoop = config_dict['use_hadoop']
        self._prod = config_dict['prod']
        self._script_list = config_dict['script_list']
        self._data_list = config_dict['data_list']
        self._generator_list = config_dict['generator_list']
        self._run_list = config_dict['run_list']
        self._ftype = config_dict['ftype']
        if config_dict.has_key('microtree'):
            self._microtree = config_dict['microtree']
        else:
            self._microtree = HIGHLAND_SCRIPTS[0]
        if config_dict.has_key('paramfile'):
            self._paramfile = config_dict['paramfile']
        else:
            self._paramfile = ''
        if config_dict.has_key('weights_store'):
            self._weights_store = config_dict['weights_store']
        else:
            self._weights_store = directories.WEIGHTS_HIGHLAND2_STORE

        self._cmd_builder = CmdBuilder()

        # flags to check whether to overwrite outputs and suppress queries
        self._overwrite = False
        self._suppress = False


    def loop(self):
        """ Where all the dirty work of putting cmds together happens
        """
        (options, args) = self.parser.parse_known_args()
        self._suppress = self._overwrite = options.force

        base_dir = (directories.HADOOP_BASE if self._use_hadoop
                    else directories.NFS_BASE)

        for a_script in self._script_list:
            for a_type in self._data_list:
                for a_gen in self._generator_list:
                    for a_run in self._run_list:
                        job_config = dict(base = base_dir,
                                          prod = self._prod,
                                          script = a_script,
                                          dtype = a_type,
                                          gen = a_gen,
                                          run = a_run,
                                          ftype = self._ftype,
                                          paramfile = self._paramfile,
                                          microtree = self._microtree,
                                          weights_store = self._weights_store)

                        data_dir_template = self._cmd_builder.data_dir_template
                        data_dir = data_dir_template.format(**job_config)

                        if not os.path.exists(data_dir):
                            print data_dir, 'does not exist!'
                        else:
                            subdirs = natsorted(os.walk(data_dir).next()[1])
                            # check if we want to process full sample or not
                            if options.dev:
                                subdirs = subdirs[0:min(len(subdirs), 4)]

                            for a_subdir in subdirs:
                                job_config['subdir'] = a_subdir
                                job_config['inp'] = os.path.join(data_dir,
                                                                 a_subdir)

                                self._append(job_config)


    def _outpath(self, cmd):
        cmd_list = shlex.split(cmd)

        # return last part of cmd which should be the dest
        # to store path for CUStageOut
        return cmd_list[-1]


    def _move_to_trash(self, source):
        """ Move source from 'store/' to '.trash/<date>/' whilst keeping
        the directory structure of source following 'store/'.

        source may be a file or a directory
        """
        date = time.strftime('%Y%m%d')
        # remove trailing slash to get parent directory
        src_dir = os.path.dirname(source.rstrip('/'))
        basepath = src_dir.replace(directories.STORE, '').strip('/')
        dest_dir = os.path.join(directories.TRASH_STORE, date, basepath)

        # now make an empty destination directory
        utils.make_dirs_if_needed(dest_dir)

        shutil.move(source, dest_dir)


    def _check_outpath(self, outpath):
        """ returns True if the outpath is clear, false otherwise
        """
        nonempty_dir = os.path.isdir(outpath) and os.listdir(outpath)
        if os.path.isfile(outpath) or nonempty_dir:
            if self._suppress:
                if self._overwrite:
                    self._move_to_trash(outpath)
            else:
                self._overwrite = utils.prompt_yes_no('Output exists: '
                                                      ''+outpath+'\nrm? ')
                self._suppress = utils.prompt_yes_no('Suppress future'
                                                     ' queries? ')
                if self._overwrite:
                    self._move_to_trash(outpath)

            return self._overwrite
        else:
            return True


    def _append(self, job):
        cmd = self._cmd_builder.build(job)

        outpath = self._outpath(cmd)
        if self._check_outpath(outpath):
            outdir, outfile = os.path.split(outpath)
            utils.make_dirs_if_needed(outdir)
            if not outfile:
                # for create_weigths there is no outfile only an outdir
                # modify the outdir path and use as job name
                key = outdir.replace(directories.STORE,
                                     '').replace('/', '_')
            else:
                key = os.path.splitext(outfile)[0]

            self.prog_dict[key] = cmd
