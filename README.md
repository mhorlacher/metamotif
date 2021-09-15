## Metamotif: Consensus Sequence Motifs (Inspired by DeepRiPe<sup>1</sup>)

### Installation
```
pip install git+https://github.com/mhorlacher/metamotif.git
```

### Example Usage

```
import metamotif

# get example data
sequence_importances = metamotif.example_data.QKI

# extract consensus motifs via embedding and clustering
meta_motifs = metamotif.extract_meta_motifs(sequence_importances)

# visualize up to 10 meta motifs
for m in meta_motifs[:10]:
    metamotif.plot_motif(m)
```

*A more detailed example can be found in the [example.ipynb](https://github.com/mhorlacher/metamotif/blob/main/setup.py/example.ipynb) notebook.*

---
### References
[1] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7050519/