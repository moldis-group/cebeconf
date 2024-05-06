import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

INSTALL_REQUIRES = [ 'pandas', 'numpy' ]

setup(
        name='cebeconf',
        version='0.1',
        packages=find_packages(),
        package_data={'cebeconf': ['data/*']},
        include_package_data=True,
        install_requires=INSTALL_REQUIRES,
)
