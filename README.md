## Getting Started
```sh
# Please download "k8" executables from github.com/attractivechaos/k8
# Generate segments
hickit.js sam2seg -v SNPs.txt sample.sam > sample.seg
# Convert .seg to .pairs (without markdup; with phase)
hickit -f -D sample.seg > sample.pairs
# Identify TADs
hickit sample.seg > sample-TAD.pairs
# Impute (without masking pairs contained in TADs)
hickit -p sample.seg > sample-im.pairs
# Plot imputed pairs
hickit -o sample.png sample-im.pairs
```
