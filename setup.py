from setuptools import setup

setup(name='structurepipeline',
      version='0.0.1',
      description='ChEMBL Structure Pipeline',
      url='https://www.ebi.ac.uk/chembl/',
      author='Greg Landrum',
      author_email='greg.landrum@t5informatics.com',
      license='MIT',
      packages=['structurepipeline'],
      package_data={'structurepipeline': ['data/*', 'data/exclude_mols/*']},
      zip_safe=False)
