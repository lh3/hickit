## Getting Started
```sh
# Download precompiled binaries for Linux
curl -L https:// | tar -jxf -
cd hickit-0.1_x64-linux

# Map Dip-C reads and extract contacts (skip if you use your own pipeline)
./pre-dip-c read1.fq.gz read2.fq.gz | bwa mem -p hs37d5.fa - | gzip > aln.sam.gz
./k8 hickit.js sam2seg -v phased_SNP.tsv aln.sam.gz | gzip > contacts.seg.gz
./hickit pair --out-phase -Dd0 contacts.seg.gz | bgzip > contacts.pairs.gz  # optional

# Impute phases (the input file can also be contacts.seg.gz)
./hickit pair -p contacts.pairs.gz | bgzip > impute.pairs.gz
./hickit pair -v.1 contacts.pairs.gz > contacts.val  # estimate phasing accuracy by holdout
# Infer 3D structure (`hickit bin -g` will get called multiple times on the input)
./fdg-multi.pl impute.pairs.gz | sh

# 2D contact map in PNG (bin size determined by the image width)
./hickit image2d -w 800 -o impute.png impute.pairs.gz
# Compute CpG density (optional; easy with your own scripts, too)
./hickit.js gfeat -r hs37d5.fa.gz impute.3dg.gz | gzip > impute.cpg.3dg.gz
# Visualize 3D structure (requiring a graphical card)
./hickit-gl view3d impute.cpg.3dg.gz
```

## Introduction

Hickit is a set of tools initially developed to process diploid single-cell
Hi-C data. It extracts contact pairs from read alignment, identifies phases of
contacts overlapping with SNPs of known phases, imputes missing phases, infers
the 3D structure of a single cell and visualizes the structure. Part of the
hickit functionality also works with bulk Hi-C data. In particular, hickit
implements a fast (untested) binning-free TAD calling algorithm and an
efficient neighboring contacts counter which can be adapted to ultrafast loop
calling.

## Installation

Hickit depends on [zlib][zlib]. The command-line tools can be compiled by
typing `make` in the source code directory. The 3D viewer (i.e. the `view3d`
command) further requires OpenGL and GLUT and can be compiled with `make gl=1`.

## Users' Guide

The "Getting Started" section above presents a brief walkthrough of the hickit
functionality. The following gives more details and explanations.

### Extracting contact pairs

[zlib]: http://zlib.net
