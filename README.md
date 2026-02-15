# Installation

```bash
$ pip install bionmf
```

# Usage

```python
from bionmf import BioNMF       # Import model

model = BioNMF()
nmf = model.fit(                # NMFInfo object, see below for properties
    adata,                      # AnnData object
    rank_range = range(2, 3),
    n_runs = 3,
    cutoff = 0.95
)

adata = model.get_adata()       # Return AnnData object
                                # Stores cell-factors under .obsm['X_nmf']
                                # Stores gene-factors under .varm['factors']

genes = model.program_genes()   # Get differential expressed genes per program
model.plot_cophcorr("cc.html")  # Cophenetic correlation for each tested rank
```

# Documentation

```python
BioNMF(
    random_state,               # Initial random state
    **kwargs                    # Rest of arguments are passed to sklearn.decomposition.NMF
).fit(
    adata,                      # Accepts AnnData and DataFrame (must be cells as rows)
    rank_range,                 # Range of rank values to test
    n_runs,                     # Number of runs (random intializations) for each rank
    cutoff                      # Cophenetic correlation cutoff
)                               # -> returns NMFInfo object

NMFInfo(                        # Properties accessed as nmf.rank, nmf.W, etc.
    rank,                       # Chosen (best) rank
    W,                          # (genes by factors) matrix -> gene programs
    H,                          # (cells by factors) matrix -> program assignment
    connectivity_mat,           # (cells by cells) matrix -> program connectivity
    reconstruction_err,         # RMSE of reconstruction
    cophcorr                    # Cophenetic correlation
)
```

# References

Brunet JP, Tamayo P, Golub TR, Mesirov JP. Metagenes and molecular pattern discovery using matrix factorization. Proc Natl Acad Sci U S A. 2004 Mar 23;101(12):4164-9. doi: 10.1073/pnas.0308531101. Epub 2004 Mar 11. PMID: 15016911; PMCID: PMC384712.
