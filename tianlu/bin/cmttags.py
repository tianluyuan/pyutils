#!/usr/bin/env python
"""
A script to create tags for CMT managed packages.

Call from within cmt/ directory
"""
import subprocess
import sys
import os
from optparse import OptionParser


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


# Ignore large external packages for now
IGNORES = ['CMT', 'EXTERN', 'GSL', 'MYSQL', 'GEANT', 'CLHEP']

# Extensions for finding src files, must satisfy unix wildcard rules
EXTENSIONS = {'cpp': ('*.[hc]', '*.[hc]xx', '*.[hc]pp', '*.cc', '*.hh'),
              'python':('*.py'),
              'java':('*.java')}

# Ignore these files and dirs, key specifies argument to find
# (e.g. '-iname')
PRUNE = {'iname':['*_Dict.[hc]*', '*linkdef.h']}


def check_dir():
    """ Are we inside cmt/
    """
    if os.path.basename(os.getcwd()) != 'cmt':
        sys.exit('Not inside cmt directory!')


def check_requirements():
    """ Ensure that requirements file exists in cmt dir
    """
    if not os.path.isfile('requirements'):
        sys.exit('No requirements file!')


def init_use_dict():
    """Returns the initial use_dict which contains the current (cwd)
    package and its path.  'cmt show uses' does not include the
    package itself.
    """
    # Must call os.path.dirname because the cwd should be inside a cmt
    # directory
    return {'this':os.path.dirname(os.getcwd())}


def parse_uses():
    """ Returns a dict of used packages and their root dir paths.
    e.g. {ROOT:/path/to/cmt/installed/ROOT/vXrY}
    """
    check_dir()
    check_requirements()
    proc = subprocess.Popen(['cmt', 'show', 'uses'],
                            stdout=subprocess.PIPE)

    use_dict = init_use_dict()
    for line in iter(proc.stdout.readline, ''):
        tokens = line.split()
        # ignore lines that start with '#'
        if line[0] != '#' and tokens[1] not in IGNORES:
            basepath = tokens[-1].strip('()')
            # highland and psyche do not strictly follow CMT path
            # organization. They have subpackages within a master, so
            # we need to take that into account
            relpath_list = [master for master in tokens[3:-1]]
            relpath_list.extend([tokens[1], tokens[2]])
            use_dict[tokens[1]] = os.path.join(basepath, *relpath_list)

    return use_dict


def get_exts(opts):
    if opts.python:
        return EXTENSIONS['python']
    elif opts.java:
        return EXTENSIONS['java']
    else:
        return EXTENSIONS['cpp']


def build_find_args(exts):
    """ ext is a list of file extensions corresponding to the files we want
    to search. This will return a list of arguments that can be passed to `find`
    """
    find_args = []
    for a_ext in exts:
        # -o for "or"
        find_args.extend(['-o', '-iname'])
        find_args.append('{0}'.format(a_ext))

    # replace first '-o' with '( for grouping matches
    find_args[0] = '('
    # append parens for grouping negation
    find_args.extend([')', '('])

    # Add prune files
    for match_type in PRUNE:
        for aprune in PRUNE[match_type]:
            find_args.append('-not')
            find_args.append('-'+match_type)
            find_args.append('{0}'.format(aprune))

    find_args.append(')')
    return find_args


def build_find_cmd(opts, paths):
    """ Builds teh cmd file using ctags.  Returns cmd based on the following

    template: 'find {0} -type f {1} | etags -'
    """
    find_args = build_find_args(get_exts(opts))

    return ['find']+paths+['-type', 'f']+find_args


def build_tags_cmd():
    return ['etags', '-']


def main():
    """ Uses ctags to generate TAGS file in cmt directory based on cmt show uses
    """
    parser = OptionParser()

    parser.add_option('--cpp',
                      dest='cpp',
                      action='store_true',
                      default=False,
                      help='tag only c/cpp files (default)')
    parser.add_option('--python',
                      dest='python',
                      action='store_true',
                      default=False,
                      help='tag only python files')
    parser.add_option('--java',
                      dest='java',
                      action='store_true',
                      default=False,
                      help='tag only java files')
    parser.add_option('-n',
                      dest='dry_run',
                      action='store_true',
                      default=False,
                      help='dry run')

    (opts, args) = parser.parse_args()

    # get the cmt show uses dictionary of programs and paths
    use_dict = parse_uses()

    # build the commands
    find_cmd = build_find_cmd(opts, list(use_dict.itervalues()))
    tags_cmd = build_tags_cmd()

    print 'Creating TAGS file based on dependencies:'
    print use_dict

    if not opts.dry_run:
        find_proc = subprocess.Popen(find_cmd, stdout=subprocess.PIPE)
        tags_proc = subprocess.Popen(tags_cmd, stdin=find_proc.stdout)

        tags_proc.communicate()


if __name__ == '__main__':
    main()
