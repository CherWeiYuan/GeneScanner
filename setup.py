from setuptools import setup, find_packages

LONG_DESCRIPTION = \
'''This program reads the output of GeneScan in csv format, remove
peaks with small or relatively small areas, and calculates, for 
each sample, the percentage of the total area that each peak covers.

The hard filtering procedure here has low threshold (in ambiguous
situations, many peaks will not be removed) and the decision to keep
or remove peaks is left to the user.'''


setup(name='genescanner',
      author='CHER_WEI_YUAN',
      version="1.0.0",
      author_email='E0031403@U.NUS.EDU',
      packages=find_packages(include=['genescanner']),
      entry_points={'console_scripts': ['genescanner=genescanner.genescanner:main']},
      url='https://github.com/CherWeiYuan/GeneScanner/',
      license='MIT LICENSE',
      description=('Bioinformatics commandline tool'),
      long_description=(LONG_DESCRIPTION),
      install_requires=['pandas==1.3.5',
                        'numpy==1.22.0',
                        'seaborn==0.11.2',
                        'matplotlib==3.5.1'])

