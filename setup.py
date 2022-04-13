from setuptools import setup

setup(
    name="chembl_structure_pipeline",
    version="1.0.1",
    description="ChEMBL Structure Pipeline",
    url="https://www.ebi.ac.uk/chembl/",
    author="Greg Landrum",
    author_email="greg.landrum@t5informatics.com",
    license="MIT",
    packages=["chembl_structure_pipeline"],
    package_data={"chembl_structure_pipeline": ["data/*"]},
    zip_safe=False,
    install_requires=["setuptools>=46.4.0"],
    extras_require={
        "rdkit": ["rdkit-pypi>=2019.03.1,<=2021.3.5.1"],
        "dev": ["bump2version<=1.0.0", "nose<=1.3.7"]
    },
)
