from setuptools import setup
import chembl_structure_pipeline

setup(
    name="chembl_structure_pipeline",
    version=chembl_structure_pipeline.__version__,
    description="ChEMBL Structure Pipeline",
    url="https://www.ebi.ac.uk/chembl/",
    author="Greg Landrum",
    author_email="greg.landrum@t5informatics.com",
    license="MIT",
    packages=["chembl_structure_pipeline"],
    package_data={"chembl_structure_pipeline": ["data/*"]},
    install_requires=['rdkit==2022.09.01'],
    zip_safe=False,
)
