## Getting Started
```sh
# Please download "k8" executables from github.com/attractivechaos/k8
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
