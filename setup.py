from setuptools import setup, find_packages

setup(name='pyutils',
      version='0.0.1',
      description='A package with useful python tools',
      url='https://github.com/tianluyuan/pyutils',
      author='Tianlu Yuan',
      author_email='tyuan@icecube.wisc.edu',
      license='MIT',
      packages=find_packages('./'),
      zip_safe=False)
