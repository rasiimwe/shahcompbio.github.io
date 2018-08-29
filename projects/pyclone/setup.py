from setuptools import find_packages, setup

setup(
    name='PyClone',
    version='0.13.1',
    description='Python tools for analysing clonal evolution using NGS data.',
    author='Andrew Roth',
    author_email='andrewjlroth@gmail.com',
    url='http://compbio.bccrc.ca',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'PyClone = pyclone.cli:main',
        ]
    }
)
