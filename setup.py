from glob import glob
from setuptools import setup, find_packages
import ngs_doit
# Run setuptools setup
setup(
    name = ngs_doit.__projectname__,
    version = ngs_doit.__version__,
    packages = find_packages(),
    scripts = glob('bin/*'),
    entry_points = {
        'console_scripts': [
            'ngs_doit_run = ngs_doit.wrapper:main'
          ]
        } )

