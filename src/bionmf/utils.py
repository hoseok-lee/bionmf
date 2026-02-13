from multiprocess import Pool
from tqdm import tqdm

import numpy as np
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import pandas as pd


def connectivity_matrix(X):

    return sum(
        # x(i) is True if belong to cluster i
        # Matrix multiplication computes to 1 if (i, j) is both True
        np.matmul(x.T, x)

        for cluster in range(X.max() + 1)
        # For each cluster membership
        # Skip if there are no elements of a certain cluster
        if (x := (X == cluster).reshape(1, -1)).any()
    )


def sparse_rmse(A, B) -> float:

    return np.sqrt(
        (A - B) \
            .power(2) \
            .mean()
    )


def log2_deg(
    adata: ad.AnnData,
    groupby: str | list,
    n_jobs: int = 1
):

    from scipy.stats import false_discovery_control
    from scanpy.tools._rank_genes_groups import _RankGenes
    from scanpy.pp import log1p, scale

    # log2 + z-scale
    log1p(adata, base = 2)
    scale(adata)
    sp_mat = adata.X

    # Passed raw data as annotations
    if not isinstance(groupby, str):
        adata.obs['groupby'] = pd.Series(
            groupby,
            index = adata.obs_names,
            dtype = "category"
        )
        groupby = 'groupby'

    # Scanpy's Wilcoxon is much faster
    # But we compute our own log foldchange
    rgg = _RankGenes(
        adata,
        groups = "all",
        groupby = groupby
    )

    # Log Foldchange
    def _proc(mask_obs):

        mat1 = g_X[mask_obs]
        mat2 = g_X[~mask_obs]

        return list(
            (
                mat1.sum(axis = 0) / mat1.shape[0] -
                mat2.sum(axis = 0) / mat2.shape[0]
            ).flat
        )


    def define_global(var):
        global g_X
        g_X = var

    # Use multiprocessing
    if n_jobs > 1:
        with Pool(
            processes = n_jobs,
            initializer = define_global,
            initargs = (sp_mat, )
        ) as pool:

            foldchange = list(
                tqdm(
                    pool.map(
                        _proc,
                        rgg.groups_masks_obs
                    ),
                    total = len(rgg.groups_masks_obs),
                    desc = "Computing log foldchanges",
                    leave = False
                )
            )

    else:
        define_global(sp_mat)
        foldchange = list(
            tqdm(
                map(
                    _proc,
                    rgg.groups_masks_obs
                ),
                total = len(rgg.groups_masks_obs),
                desc = "Computing log foldchanges",
                leave = False
            )
        )

    wilcoxon = list(
        tqdm(
            rgg.wilcoxon(tie_correct = True),
            total = len(rgg.groups_masks_obs),
            desc = "Statistical test",
            leave = False
        )
    )

    for fc, w in zip(foldchange, wilcoxon):

        df = pd.DataFrame(
            {
                'foldchange': fc,
                'pvals': w[2]
            },
            index = adata.var_names.rename("genes")
        )

        df['pvals_adj'] = false_discovery_control(
            df['pvals'],
            method = 'bh'
        )

        yield rgg.groups_order[w[0]], df
