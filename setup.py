from setuptools import setup

setup(
    name="chembl_structure_pipeline",
    version="1.0.0",
    description="ChEMBL Structure Pipeline",
    url="https://www.ebi.ac.uk/chembl/",
    author="Greg Landrum",
    author_email="greg.landrum@t5informatics.com",
    license="MIT",
    packages=["chembl_structure_pipeline"],
    package_data={"chembl_structure_pipeline": ["data/*"]},
    zip_safe=False,
    install_requires=[
          "rdkit-pypi>=2019.03.1,<=2022.3.1"
    ]
)
