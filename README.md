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
# 2D contact map in PNG (bin size determined by the image width)
./hickit image2d -w 800 -o impute.png impute.pairs.gz

# Infer 3D structure (which calls `hickit bin -g` multiple times on the input)
./fdg-multi.pl impute.pairs.gz | sh
# Compute CpG density (optional)
./hickit.js gfeat -r hs37d5.fa.gz impute.3dg.gz | gzip > impute.cpg.3dg.gz
# Visualize 3D structure (requiring a graphical card)
./hickit-gl view3d impute.cpg.3dg.gz
```
