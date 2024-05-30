

### 6. Quality Control

It might be worth considering automatic thresholding via MAD (median absolute deviations). The MAD is given by 
$$MAD = median(X_i - median(X))$$
with $X_i$ being the respective QC metric of an observation and describes a robust statistic of the variability of the metric. we mark cells as outliers if they differ by 5 MADs which is a relatively permissive filtering strategy.

```
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
```

It is important to note that mitochondrial counts are annotated either with the prefix “mt-” or “MT-’’ depending on the species considered in the dataset. As mentioned, the dataset used in this notebook is human bone marrow, so mitochondrial counts are annotated with the prefix “MT-”. For mouse datasets, the prefix is usually lower case, so “mt-“.

For droplet based single cell RNA-seq experiments, a certain amount of background mRNAs is present in the dilution that gets distributed into the droplets with cells and sequenced along with them. The net effect of this is to produce a background contamination that represents expression not from the cell contained within a droplet, but the solution that contained the cells. Cell-free mRNA molecules, also known as ambient RNA, can confound the number of observed counts and can be seen as background contamination.

SoupX calculates the profile of the soup. It estimates the ambient mRNA expression profile from empty droplets as given by the unfiltered Cellranger matrix. Next, SoupX estimates the cell specific contamination fraction. Lastly, it corrects the expression matrix according to the ambient mRNA expression profile and the estimated contamination. The output of SoupX is a modified counts matrix.

We additionally filter out genes that are not detected in at least 20 cells as these are not informative.

A doublet is called homotypic if it is formed by the same cell type (but from different individuals) and heterotypic otherwise. Homotypic doublets are not necessarily identifiable from count matrices and are often considered innocuous as they can be identified with cell hashing or SNPs. Hence, their identification is not the main goal of the doublet detection methods.

Doublets formed from different cell types or states are called heterotypic. Their identification is crucial as they are most likely misclassified and can lead to distorted downstream analysis steps.

![double detection](https://www.sc-best-practices.org/_images/doublet_detection.jpeg)
Common doublet detection methods generate artificial doublets by randomly subsampling pairs of cells and averaging their gene expression profile to obtain doublet counts.

Key takeaways:
1. Filtering of poor-quality cells should be based on median absolute deviations with lenient cutoffs to avoid bias against smaller subpopulations.
2. Feature-based filtering does not show benefits for downstream tasks.
3. Doublets can be efficiently detected with tools like scDblFinder.
4. Doublet detection methods should not be run on aggregated scRNA-seq data representing multiple batches.
