# microbiome_structure_v3
this is a reiteration of microbiome_structure_v2 with cleaner folder tree and fixed intralayer edge weights.

This repository is accompanying the paper "Scale-dependent signatures of microbial co-occurrence revealed via multilayer network analysis". 
The paper is still in review and made available for reviewers examination. It is not yet intended for public access.

# Cow genetic analysis

when working with FWASS (see methods):

Code for creating the genetic similarity matrix: 
\texttt{python3 compute\_similarity.py -i INPUT\_FILE.vcf --weighted --ns -o OUTPUT\_FILE}

the “--weighted” flag calculates an allele frequency-weighted metric 
the “—ns” flag converts the output into a format that is compatible with the input format required for Netstruct Hierarchy (.txt file)

Source:  https://github.com/Greenbaum-Lab/FWASS.git

See https://github.com/Greenbaum-Lab/FWASS for further details.

