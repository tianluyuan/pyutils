from setuptools import setup, find_packages

setup(name='pyutils_hep',
      version='0.1.3',
      description='A package with useful python tools',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      url='https://github.com/tianluyuan/pyutils',
      author='Tianlu Yuan',
      author_email='tyuan@icecube.wisc.edu',
      license='MIT',
      packages=find_packages('./'),
      install_requires=['numpy>=1.17',
                        'scipy'],
      extras_require={
          'plotting':  ['matplotlib']},
      )
