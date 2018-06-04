## Getting Started
```sh
# Download precompiled binaries for Linux
curl -L https:// | tar -jxf -
cd hickit-0.1_x64-linux
# Map Dip-C reads and extract contacts (skip if you used your own pipeline)
./pre-dip-c read1.fq.gz read2.fq.gz | bwa mem -p human.fa - | gzip > aln.sam.gz
./k8 hickit.js sam2seg -v phased_SNP.tsv aln.sam.gz | gzip > contacts.seg.gz
# 
# Generate segments
hickit.js sam2seg -v SNPs.txt sample.sam > sample.seg
# Convert .seg to .pairs (without markdup; with phase)
hickit --out-phase -D sample.seg > sample.pairs
# Identify TADs (both .seg and .pairs work)
hickit -t sample.seg > sample-TAD.pairs
hickit -t sample.pairs > sample-TAD.pairs
# Impute (without masking pairs contained in TADs)
hickit -p sample.seg > sample-im.pairs
# Plot imputed pairs
hickit -I sample.png sample-im.pairs
```
