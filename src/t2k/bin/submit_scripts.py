#!/usr/bin/env python
#
# submit shell script to cluster
import os
from tianlu.submit import submit
from tianlu import utils


__author__ = 'Tianlu Yuan'
__email__ = 'tianlu.yuan [at] colorado.edu'


class submitScripts(submit):
    def loop(self):
        (options, args) = self.parser.parse_args()
        sh_file_list = args

        for sh_file in sh_file_list:
            key = os.path.splitext(os.path.basename(sh_file))[0]

            self.prog_dict[key] = utils.abspath(sh_file)



if __name__ == '__main__':
    submitScripts().run()

