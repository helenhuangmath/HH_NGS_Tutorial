# 15 · Spatial transcriptomics: integrating sections with Scanpy + Harmony

[← scATAC-seq](14_scatacseq.md) · [Index](../README.md)

Several tissue sections or regions of interest (ROIs) are clustered **jointly**, so that a cluster means the same tissue domain in every section. Batch effects between sections are corrected with Harmony before building the neighbour graph.

```text
per-section AnnData (.h5ad, with .obsm["spatial"])
  └─► concat ─► normalise / log ─► HVGs ─► PCA ─► Harmony ─► neighbours (X_pca_harmony)
        └─► Leiden clusters ─► UMAP ─► clusters on tissue ─► markers / composition / neighbourhoods
```

---

## 15.1 Integrate sections

```python
import scanpy as sc, scanpy.external as sce, pandas as pd
from datetime import date

sc.set_figure_params(dpi_save=400, facecolor="white", figsize=(4, 4), fontsize=7)
tag = date.today().strftime("%Y%m%d")

samples = ["WT_S1_roi1", "WT_S1_roi2", "WT_S2_roi1", "WT_S2_roi2",
           "KO_S1_roi1", "KO_S1_roi2", "KO_S2_roi1", "KO_S2_roi2"]

adatas = {}
for s in samples:
    a = sc.read_h5ad(f"data/{s}/adata_{s}.h5ad")           # load into memory (not backed)
    a.var_names_make_unique()
    adatas[s] = a

adata = sc.concat(adatas, label="sample", uns_merge="unique", index_unique="-")
adata.obs["sample"]   = pd.Categorical(adata.obs["sample"], categories=samples)
adata.obs["genotype"] = adata.obs["sample"].str.split("_").str[0]
adata.layers["counts"] = adata.X.copy()

# QC
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.filter_cells(adata, min_counts=50)
sc.pp.filter_genes(adata, min_cells=10)

# normalise → log → PCA
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pp.pca(adata, n_comps=50, mask_var="highly_variable")

# Harmony corrects the PCs → stored in .obsm["X_pca_harmony"]
sce.pp.harmony_integrate(adata, key="sample")

# build the graph on the corrected embedding
sc.pp.neighbors(adata, use_rep="X_pca_harmony")
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters", resolution=1.0)

adata.write("ADATA_integrated.h5ad")
```

> **Common pitfalls**
> 1. `sc.pp.neighbors()` uses `X_pca` by default. Without `use_rep="X_pca_harmony"`, the UMAP and clusters are **not** integrated.
> 2. Normalise and log-transform once, before PCA and Harmony.
> 3. Objects opened with `backed="r"` cannot be concatenated or modified in memory; read them normally.
> 4. Use the same `obs` key (here `sample`) for concatenation, Harmony and plotting.
>
> `sc.pp.normalize_per_cell` is deprecated; use `sc.pp.normalize_total`. In scanpy < 1.10, use `use_highly_variable=True` instead of `mask_var=`.

---

## 15.2 Visualise

```python
import matplotlib.pyplot as plt

# before vs after integration (run PCA-based UMAP once without Harmony for the "before" view)
sc.pl.umap(adata, color=["sample", "clusters"], wspace=0.4, save=f"_afterHarmony_{tag}.pdf", show=False)

# clusters on tissue coordinates, one panel per section
fig, axes = plt.subplots(2, 4, figsize=(7.2, 3.6))
for ax, s in zip(axes.ravel(), samples):
    sc.pl.embedding(adata[adata.obs["sample"] == s], basis="spatial", color="clusters",
                    title=s, ax=ax, show=False, legend_loc=None, s=8)
    ax.set_aspect("equal"); ax.axis("off")
fig.savefig(f"figures/spatial_clusters_{tag}.pdf", bbox_inches="tight")

# a gene on tissue
sc.pl.embedding(adata[adata.obs["sample"] == "KO_S1_roi1"], basis="spatial", color="GeneX",
                cmap="magma_r", vmax="p99", save=f"_GeneX_KO_S1_{tag}.pdf", show=False)

# cluster composition per section
comp = pd.crosstab(adata.obs["sample"], adata.obs["clusters"], normalize="index")
comp.plot(kind="barh", stacked=True, width=0.7, figsize=(3.5, 2.2), legend=False)
plt.xlabel("Fraction of spots"); plt.tight_layout(); plt.savefig(f"figures/composition_{tag}.pdf")
```

Markers per cluster and genotype comparisons:

```python
sc.tl.rank_genes_groups(adata, groupby="clusters", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, save=f"_markers_{tag}.pdf", show=False)
```

For KO vs WT within a cluster, use pseudobulk counts per section (`adata.layers["counts"]` summed by `sample` × `clusters`) with DESeq2 (chapter 5); spots from one section are not independent replicates.

### Spatial statistics (Squidpy)

```python
import squidpy as sq
sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial", library_key="sample")
sq.gr.nhood_enrichment(adata, cluster_key="clusters")     # which clusters are neighbours
sq.pl.nhood_enrichment(adata, cluster_key="clusters", cmap="RdBu_r")
sq.gr.spatial_autocorr(adata, mode="moran")                # spatially variable genes (Moran's I)
```

---

## Example output

![Fig. 15](figures/fig15_spatial.png)

**Fig. 15 | Integration of spatial transcriptomics sections.** **a**, UMAP before integration coloured by genotype; sections separate by batch. **b**, UMAP after Harmony coloured by Leiden cluster. **c**, Cluster composition per section. **d**,**e**, Clusters mapped back to tissue coordinates in a WT and a KO section; a domain (cluster 4) expands in the KO section. **f**, Expression of a cluster-4 marker on the KO section. Scale bar, 500 µm. *Simulated data for illustration.*

---

## References and resources

- Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* **19**, 15 (2018). — <https://scanpy.readthedocs.io>
- Korsunsky, I. *et al.* Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat. Methods* **16**, 1289–1296 (2019). — `harmonypy`: <https://github.com/slowkow/harmonypy>
- Palla, G. *et al.* Squidpy: a scalable framework for spatial omics analysis. *Nat. Methods* **19**, 171–178 (2022). — <https://squidpy.readthedocs.io>
- Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. *Sci. Rep.* **9**, 5233 (2019).
- Virshup, I. *et al.* The scverse project provides a computational ecosystem for single-cell omics data analysis. *Nat. Biotechnol.* **41**, 604–606 (2023). — <https://scverse.org>
- Single-cell best practices, spatial chapters — <https://www.sc-best-practices.org>

---

[← scATAC-seq](14_scatacseq.md) · [Index](../README.md)
