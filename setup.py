from setuptools import setup

__version__ = '0.0.1'

setup(name='structurepipeline',
      version=__version__,
      description='ChEMBL Structure Pipeline',
      url='https://www.ebi.ac.uk/chembl/',
      author='Greg Landrum',
      author_email='greg.landrum@t5informatics.com',
      license='MIT',
      packages=['structurepipeline'],
      package_data={'structurepipeline': ['data/*']},
      zip_safe=False)
