from setuptools import setup, find_packages
  
# reading long description from file
with open('DESCRIPTION.txt') as file:
    long_description = file.read()
  

# some more details
CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'Environment :: MacOS X',
    'Framework :: IDLE',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Natural Language :: English',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.11',
    'Topic :: Scientific/Engineering :: Bio-Informatics'
    'Topic :: Scientific/Engineering :: Information Analysis'
    ]
  
# calling the setup function 
setup(name='Chimera_Buster',
      version='1.3.0',
      description='eliminates chimeric reads by comparing the UMI sequences and finding any matches in the 5 prime or 3 prime UMIs and keeping the sequence that has the highest prevalence',
      long_description=long_description,
      url='https://github.com/JessicaA2019/Chimera_Buster',
      author='Jessica Lauren ALbert',
      author_email='jessica.albert2001@gmail.com',
      license='MIT',
      packages = find_packages(),
      entry_points = {'console_scripts': ['Chimera_Buster = Chimera_buster.CLI_chimera_buster:main']},
      classifiers=CLASSIFIERS,
      keywords='HIV UMIs PCR Chimeras ONTsequencing',
      include_package_data = True
      )

