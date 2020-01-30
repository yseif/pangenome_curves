import os
from setuptools import setup, Command, find_packages


class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')
setup(
    name='pangenome_curves',
    version='0.0.1',
    description='A package to draw pan genome curves.',
    license='MIT',
    packages=find_packages(),
    long_description=open('README.md').read(),
    url='https://github.com/yseif/pangenome_curves/',
    author='Yara Seif',
    cmdclass={
        'clean': CleanCommand,
    }
)
