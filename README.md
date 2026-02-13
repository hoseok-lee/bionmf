# Installation

```bash
$ pip install bionmf
```

# Usage

```python
from bionmf import BioNMF       # Import model

model = BioNMF()
nmf = model.fit(                # NMFInfo object
    adata,                      # AnnData object
    rank_range = range(2, 3),
    n_runs = 3,
    cutoff = 0.95
)

genes = model.program_genes()   # Get differential expressed genes per program

model.plot_cophcorr()           # Cophenetic correlation for each tested rank
model.plot_heatmap()            # Heaetmap of gene programs
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
    cutoff,                     # Cophenetic correlation cutoff
)                               # -> returns NMFInfo object

NMFInfo(
    rank,                       # Best rank
    W,                          # (genes by factors) matrix -> gene programs
    H,                          # (cells by factors) matrix -> program assignment
    connectivity_mat,           # (cells by cells) matrix -> program connectivity
    reconstruction_err,         # RMSE of reconstruction
    cophcorr,                   # Cophenetic correlation
)
```
