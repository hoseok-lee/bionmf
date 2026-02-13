from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import squareform
from fastcluster import linkage
from collections import defaultdict
from typing import Literal
from dataclasses import dataclass
from typing import Optional, Any
from anndata import AnnData
from pathlib import Path

import scipy.sparse as sp
import seaborn as sns
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

from .utils import sparse_rmse, log2_deg, connectivity_matrix


class BioNMF:

    @dataclass
    class NMFInfo:

        rank: int
        W: np.ndarray
        H: np.ndarray
        connectivity_mat: float
        reconstruction_err: float
        cophcorr: float

        @property
        def programs(self):
            return np.argmax(self.H, axis = 0)


    def __init__(
        self,
        random_state: int = 0,
        **kwargs
    ):

        # Annotations for cell and gene names
        # Ideally, this is an AnnData object
        self.data = None
        self.obs_names = None
        self.var_names = None

        self.random_state = random_state
        self.nmf_runs = defaultdict(self.NMFInfo)

        # Default NMF solver arguments
        self.nmf_kwargs = {
            'init': "random"
        } | kwargs


    def _select_rank(
        self,
        nmf_runs: dict,
        cutoff: float = 0.95
    ):

        # Rank is selected by the following two conditions
        # 1. First rank with cophenetic correlation below cutoff after two
        #    previous ranks above threshold
        # 2. Closest correlation to cutoff adjacent to crossover point
        # For example, with cutoff 0.95
        # [ 1, 0.98, 0.93, 0.97, 0.97, 0.94, 0.93, ... ]
        #                                ^ selects this rank

        # Compute the difference of cophenetic correlation from cutoff
        # Positive values above cutoff, negatives below cutoff
        diff = np.array([nmf_runs[rank].cophcorr for rank in nmf_runs]) - cutoff

        # Let's dissect this...
        arg = np.argwhere(
            np.convolve(
                # Converts the difference into a binary array such that
                # positive values (and 0) are 1, negative values are -1
                (diff >= 0) * 2 - 1,

                # Through convolution, seeks a pattern such that there is a
                # positive, positive, negative difference (above, above, below)
                # (Note that convolution will flip this before mapping)
                [ -1, 1, 1 ]
            # Convolution value will equal 3 if pattern is matched
            ) == 3
        )

        # No pattern found, just go with max rank
        if not np.any(arg):
            return max(nmf_runs)

        # Choose the first crossing point, there can be multiple
        cross = arg.flatten()[0]
        idx = cross if abs(diff[cross]) < abs(diff[cross - 1]) else cross - 1

        return list(nmf_runs.keys())[idx]


    def _extract_X(self, data):

        # AnnData object -> assume .X stores data
        if isinstance(data, AnnData):
            self.X = data.X.T
            self.adata = data

        # pd.DataFrame
        if isinstance(data, pd.DataFrame):
            self.X = data.T.values
            self.adata = AnnData(data)

        # Accepts np.matrix, sp.spmatrix
        else:
            self.X = data.T
            self.adata = AnnData(data)

        return self.X


    def _fit_check(self):

        # Check if the model was fit at all
        if len(self.nmf_runs) == 0:
            raise ValueError(
                "Please make sure to fit the model on at least one NMF rank and"
                "one run with BioNMF.fit()"
            )


    def fit(
        self,
        data: Any,
        rank_range: list = range(2, 20 + 1),
        n_runs: int = 5,
        cutoff: float = 0.95
    ):

        # Extract matrix
        self.X = self._extract_X(data)
        self.n_runs = n_runs
        self.rank_range = rank_range
        self.cutoff = cutoff


        for rank in self.rank_range:

            # Row by row <- cummulative connectivity matrix
            conns = np.zeros((self.X.shape[1], self.X.shape[1]))
            best_fit = -1

            W = np.zeros((self.X.shape[0], rank))
            H = np.zeros((rank, self.X.shape[1]))

            for run in range(self.n_runs):

                model = NMF(
                    n_components = rank,
                    random_state = self.random_state + run,
                    **self.nmf_kwargs
                )

                W_ = model.fit_transform(self.X)
                H_ = model.components_

                # Use RMSE, reconstruction_err_ from sklearn uses Frobenius
                rmse = sparse_rmse(
                    self.X,
                    sp.csr_matrix(np.matmul(W_, H_))
                )

                # Cell state connectivity
                cell_states = np.argmax(H_, axis = 0)
                conns += connectivity_matrix(cell_states)

                # Although it is ugly, minimizes number of loops
                # Update best fitting NMF
                if (best_fit == -1) or (rmse < best_fit):
                    best_fit = rmse
                    W = W_
                    H = H_

            # Consensus connectivity matrix
            C = conns / self.n_runs

            # Convert distance (1 - C) to condensed matrix form
            # squareform performs inverse if square matrix --> upper-triangle
            d = squareform(1 - C)

            cophcorr, _ = cophenet(linkage(d, method = 'average'), d)

            self.nmf_runs[rank] = self.NMFInfo(
                rank = rank,
                W = W,
                # Normalize by column --> each cell has probability of state
                H = H / H.sum(axis = 0),
                connectivity_mat = C,
                reconstruction_err = best_fit,
                cophcorr = cophcorr
            )

        # Choose rank based on Brunet et al method
        self.rank = self._select_rank(self.nmf_runs, cutoff = cutoff)

        return self.nmf_runs[self.rank]


    def plot_cophcorr(
        self,
        save: Optional[str] = None
    ):

        self._fit_check()

        # Cophenetic coefficient
        plt.clf()
        plt.plot(
            self.nmf_runs.keys(),
            [
                self.nmf_runs[rank].cophcorr
                for rank in self.nmf_runs
            ],
            marker = 'o',
            color = 'black',
            linestyle = '-',
            label = "Cophenetic Coefficient"
        )

        # Cutoff point
        plt.axhline(
            y = self.cutoff,
            color = 'black',
            linestyle = '--',
            label = "Cophenetic Cutoff"
        )

        # Selected rank
        plt.axvline(
            x = self.rank,
            color = 'red',
            linestyle = '--',
            label = "Selected Rank"
        )

        plt.xlabel("Rank")
        plt.ylabel("Cophenetic Coefficient")

        # This is strictly to avoid floating point values being labeled for
        # x-ticks --> force them to be integers
        plt.xticks(list(map(int, self.nmf_runs)))
        plt.legend()
        plt.tight_layout()
        plt.show()

        if save:
            plt.savefig(save)


    def program_genes(
        self,
        pval_cutoff: Optional[float] = None,
    ):

        programs = pd.Series(
            # Best NMF from fitting
            self.nmf_runs[self.rank] \
                .programs \
                .astype(str),

            index = self.adata.obs_names,
            name = "NMF Program"
        ).sort_values()

        # Organize genes based on differential expression per state
        results = defaultdict(pd.DataFrame)

        for program, df in log2_deg(
            self.adata,
            groupby = programs
        ):

            df = df[df['foldchange'] > 0]
            if pval_cutoff: # Cutoff for significance threshold
                df = df[df['pvals_adj'] < pval_cutoff]
            results[program] = df

        return pd.concat(
            results.values(),
            axis = 0,
            keys = results,
            names = ['program', 'genes']
        )


    def plot_heatmap(
        self,
        save: Optional[str] = None
    ):

        df = self.program_genes()
        self.adata = self.adata[
            :, df \
                .reset_index() \
                .loc[lambda x: x.groupby(by = "genes")['foldchange'].idxmax()] \
                .sort_values(by = "program")['genes']
        ]

        # Columns sorted by cell state
        # The reason for doing this operation after the genes is to make the
        # check at line 258 (index equality check) more efficient
        self.adata = self.adata[programs.index, :]

        lut = dict(
            zip(
                programs.unique(),
                plt.get_cmap('Set1')(
                    np.linspace(
                        0, 1,
                        programs.nunique()
                    )
                )
            )
        )
        col_colors = programs.map(lut)

        plt.clf()
        g = sns.clustermap(
            self.adata.to_df().T,
            cmap = 'bwr',
            z_score = 0,
            center = 0,
            vmin = -2, vmax = 2,
            xticklabels = False,
            yticklabels = False,
            row_cluster = False,
            col_cluster = False,
            col_colors = col_colors,
            cbar_pos = (0.09, 0.2, 0.03, 0.4),
            figsize = (5.5, 5)
        )

        # Remove axes
        g.ax_heatmap.set_xticklabels([])
        g.ax_heatmap.set_yticklabels([])
        g.ax_heatmap.set_xlabel("")
        g.ax_heatmap.set_ylabel("")

        x0, *_ = g.cbar_pos
        g.ax_cbar.set_title(
            "Relative Expression",
            rotation = 90,
            x = x0 - 1.2,
            y = 0.15,
            fontsize = 10
        )

        g.ax_col_dendrogram.legend(
            [
                plt.Rectangle(
                    (0, 0),
                    1, 1,
                    fc = color,
                    label = f"NMF{program}"
                )
                for program, color in lut.items()
            ],
            list(map(lambda x: f"NMF{x}", lut)),
            title = 'NMF Program',
            loc = 'center',
            bbox_to_anchor = (0.5, 0.6),
            frameon = False,
            ncol = min(5, len(lut))
        )

        plt.tight_layout()
        plt.show()

        if save:
            plt.savefig(save)

        return df
