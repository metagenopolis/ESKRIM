from setuptools import find_packages, setup
from os import path

cur_directory = path.abspath(path.dirname(__file__))
with open(path.join(cur_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

setup(
    name="eskrim",
    version="1.0.0",
    license="GPLv3",
    description="ESKRIM is a reference-free tool that compares microbial richness in shotgun metagenomic samples by counting k-mers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires=">=3.6",
    url="https://forgemia.inra.fr/florian.plaza-onate/eskrim",
    keywords=["bioinformatics", "metagenomics", "microbial richness", "k-mers"],
    author="Florian Plaza Oñate, Emmanuelle Le Chatelier", 
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Software Development :: Libraries",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GPLv3 License",
        "Programming Language :: Python :: 3",
    ],
)
