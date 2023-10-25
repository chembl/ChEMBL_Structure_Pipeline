from setuptools import setup

setup(
    name="chembl_structure_pipeline",
    description="ChEMBL Structure Pipeline",
    url="https://www.ebi.ac.uk/chembl/",
    author="Greg Landrum, Eloy Félix",
    author_email="greg.landrum@t5informatics.com",
    license="MIT",
    packages=["chembl_structure_pipeline"],
    package_data={"chembl_structure_pipeline": ["data/*"]},
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    install_requires=["setuptools>=46.4.0", "rdkit>=2022.09.01"],
    zip_safe=False,
)
