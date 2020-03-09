# ESKRIM: EStimate with K-mers the RIchness in a Microbiome #

### Purpose ###

* ESKRIM is a reference-free tool that compares microbial richness in metagenomic samples by counting k-mers

### Requirements ###
* Python3
* [Jellyfish2](https://github.com/gmarcais/Jellyfish) with Python bindings

### Usage ###
```
python3 eskrim.py -h
```

Contaminants reads from the host genome should be removed before using ESKRIM.
If your FASTQ files are paired, we recommend to use only forward reads.

### Authors ###

* Florian Plaza Oñate: florian.plaza-onate@inrae.fr
* Emmanuelle Le Chatelier: emmanuelle.le-chatelier@inrae.fr

