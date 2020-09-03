from setuptools import find_packages, setup

setup(
    name="eskrim",
    version="1.0.0",
    license="GPLv3",
    description="ESKRIM is a reference-free tool that compares microbial richness in shotgun metagenomic samples by counting k-mers",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    python_requires=">=3.6",
    url="https://forgemia.inra.fr/florian.plaza-onate/eskrim",
    keywords=["bioinformatics", "metagenomics", "microbial richness", "k-mers"],
    author="Florian Plaza Oñate, Emmanuelle Le Chatelier" 
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Software Development :: Libraries",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GPLv3 License",
        "Programming Language :: Python :: 3",
    ],
)
