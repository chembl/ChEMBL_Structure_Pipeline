import os
from setuptools import setup

version_file = open('VERSION')
version = version_file.read().strip()

setup(name='chembl_structure_pipeline',
      version=version,
      description='ChEMBL Structure Pipeline',
      url='https://www.ebi.ac.uk/chembl/',
      author='Greg Landrum',
      author_email='greg.landrum@t5informatics.com',
      license='MIT',
      packages=['chembl_structure_pipeline'],
      package_data={'chembl_structure_pipeline': ['data/*']},
      zip_safe=False)
