import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

INSTALL_REQUIRES = [ 'pandas', 'numpy' ]

setup(
    name='cebeconf',
    version='1.0.1',
    packages=find_packages(),
    package_data={'cebeconf': ['data/*']},
    author='Raghunathan Ramakrishnan'
    author_email='raghu.rama.chem@gmail.com'
    include_package_data=True,
    url='https://github.com/moldis-group/cebeconf'
    license='MIT License'
    description='cebeconf: A package of machine-learning models for predicting 1s-core electron binding energies of CONF atoms in organic molecules.'
    long_desc_type="text/markdown"
    install_requires=[ 'pandas', 'numpy' ]
    include_package_data=True,
)

