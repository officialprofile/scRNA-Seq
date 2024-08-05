

## Introduction

### 2. Single-cell RNA sequencing

Bulk RNA-Seq results in cell-averaged expression profiles, which are generally easier to analyze, but also hide some of the complexity such as cell expression profile heterogeneity, which may help answer the question of interest. Some drugs or perturbations may affect only specific cell types or interactions between cell types.

There are two major approaches to transcript quantification: full-length and tag-based. Full-length protocols try to cover the whole transcript uniformly with sequencing reads, whereas tag-based protocols only capture the 5’ or 3’ ends. 

Full-length sequencing is restricted to plate-based protocols (see below) and the library preparation is comparable to bulk RNA-seq sequencing approaches. An even coverage of transcripts is not always achieved with full-length protocols and therefore specific regions across the gene body may still be biased. A major advantage of full-length protocols is that they allow for the detection of splice variants. Tag-based protocols only sequence the 3’ or 5’ ends of the transcripts. This comes at the cost of not (necessarily) covering the full gene length, making it difficult to unambiguously align reads to a transcript and distinguishing between different isoforms. However, it allows for the usage of unique molecular identifiers (UMIs), which are useful to resolve biases in the transcript amplification process. The transcript amplification process is a critical step in any RNA-seq sequencing run, to ensure that the transcripts are abundant enough for quality control and sequencing. During this process, which is typically conducted with polymerase chain reaction (PCR), copies are made from identical fragments of the original molecule. Since the copies and the original molecules are indistinguishable, determining the original number of molecules in samples becomes challenging. The usage of UMIs is a common solution to quantify the original, non-duplicated molecules. The UMIs serve as molecular barcodes and are also sometimes referred to as random barcodes. These ‘barcodes’ consist of short random nucleotide sequences that are added to every molecule in the sample as a unique tag. UMIs must be added during library generation before the amplification step.

Microfluidic device based single-cell strategies trap cells inside hydrogel droplets allowing for compartmentalisation into single-cell reaction chambers. The most widely used protocols inDrop, Drop-seq and the commercially available 10x Genomics Chromium. In microfluidic based protocols only about 10% of the transcripts of the cell are recovered. Notably, this low sequencing is sufficient for robust identification of cell types. A comparison from Zhang et al. in 2019 uncovered that inDrop and Drop-seq are outperformed by 10X Genomics with respect to bead quality, as the cell barcodes in the former two systems contained obvious mismatches. Moreover, the proportion of reads originating from valid barcodes was 75% for 10X Genomics, compared to only 25% for InDrop and 30% for Drop-seq.

Similar advantages were demonstrated for 10X Genomics regarding sensitivity. During their comparison, 10X Genomics captured about 17000 transcripts from 3000 genes on average, compared to 8000 transcripts from 2500 genes for Drop-seq and 2700 transcripts from 1250 genes for InDrop. Although 10X Genomics was shown to outperform the other protocols in various aspects, it is also about twice as expensive per cell. Limitations: Captures 3’ only and not full transcripts, because the cell barcodes and PCR handles are only added to the end of the transcript.

Plate based protocols typically separate the cells physically into microwell plates. The first step entails cell sorting by, for example, fluorescent-activated cell sorting (FACS), where cells are sorted according to specific cell surface markers; or by micro pipetting. Depending on the protocol, plate based protocols might be labor-intensive with many required pipetting steps, leading to potential technical noise and batch effects.

Single-cell profiling does not always provide an unbiased view on cell types for specific tissues or organs, such as, for example, the brain. oreover, single-cell sequencing highly relies on fresh tissue, making it difficult to make use of tissue biobanks. On the other hand, the nuclei are more resistant to mechanical force, and can be safely isolated from frozen tissue without the use of tissue dissociation enzymes.

### 3. Raw data processing
### 4. Analysis frameworks and tools
### 5. Interoperability

## Preprocessing and visualization

### 6. Quality Control

It might be worth considering automatic thresholding via **MAD (median absolute deviations)**. The MAD is given by 
$$MAD = median(X_i - median(X))$$
with $X_i$ being the respective QC metric of an observation and describes a robust statistic of the variability of the metric. we mark cells as outliers if they differ by 5 MADs which is a relatively permissive filtering strategy.

```python
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
```

It is important to note that mitochondrial counts are annotated either with the prefix **“mt-” or “MT-’’** depending on the species considered in the dataset. As mentioned, the dataset used in this notebook is human bone marrow, so mitochondrial counts are annotated with the prefix “MT-”. For mouse datasets, the prefix is usually lower case, so “mt-“.

For droplet based single cell RNA-seq experiments, a certain amount of **background mRNAs** is present in the dilution that gets distributed into the droplets with cells and sequenced along with them. The net effect of this is to produce a background contamination that represents expression not from the cell contained within a droplet, but the solution that contained the cells. Cell-free mRNA molecules, also known as ambient RNA, can confound the number of observed counts and can be seen as background contamination.

SoupX calculates the profile of the soup. It estimates the ambient mRNA expression profile from empty droplets as given by the unfiltered Cellranger matrix. Next, SoupX estimates the cell specific contamination fraction. Lastly, it corrects the expression matrix according to the ambient mRNA expression profile and the estimated contamination. The output of SoupX is a modified counts matrix.

We additionally filter out genes that are not detected in at least 20 cells as these are not informative.

A doublet is called **homotypic** if it is formed by the same cell type (but from different individuals) and **heterotypic** otherwise. Homotypic doublets are not necessarily identifiable from count matrices and are often considered innocuous as they can be identified with cell hashing or SNPs. Hence, their identification is not the main goal of the doublet detection methods. Doublets formed from different cell types or states are called heterotypic. Their identification is crucial as they are most likely misclassified and can lead to distorted downstream analysis steps.

![double detection](https://www.sc-best-practices.org/_images/doublet_detection.jpeg)
Common doublet detection methods generate artificial doublets by randomly subsampling pairs of cells and averaging their gene expression profile to obtain doublet counts.

Key takeaways:
1. Filtering of poor-quality cells should be based on median absolute deviations with lenient cutoffs to avoid bias against smaller subpopulations.
2. Feature-based filtering does not show benefits for downstream tasks.
3. Doublets can be efficiently detected with tools like scDblFinder.
4. Doublet detection methods should not be run on aggregated scRNA-seq data representing multiple batches.


### 7. Normalization

The **shifted logarithm** works beneficial for stabilizing variance for subsequent dimensionality reduction and identification of differentially expressed genes. **Scran** was extensively tested and used for batch correction tasks and **analytic Pearson residuals** are well suited for selecting biologically variable genes and identification of rare cell types.

### 8. Feature selection

Usually, the scRNA-seq experiment and resulting dataset focuses on one specific tissue and hence, only a small fraction of genes is informative and biologically variable. Traditional approaches and pipelines either compute the coefficient of variation (highly variable genes) or the average expression level (highly expressed genes) of all genes to obtain 500-2000 selected genes and use these features for their downstream analysis steps. However, these methods are highly sensitive to the normalization technique used before.

### 9. Dimensionality Reduction

Xiang et al. compared in an independent comparison the stability, accuracy and computing cost of 10 different dimensionality reduction methods [Xiang et al., 2021]. They propose to use t-distributed stochastic neighbor embedding (t-SNE) as it yielded the best overall performance. Uniform manifold approximation and projection (UMAP) showed the highest stability and separates best the original cell populations.

## Identifying cellular structure

### 10. Clustering

In scRNA-seq data analysis, we describe cellular structure in our dataset with finding cell identities that relate to known cell states or cell cycle stages. This process is usually called cell identity annotation. For this purpose, we structure cells into clusters to infer the identity of similar cells. The Leiden algorithm is as an improved version of the Louvain algorithm which outperformed other clustering methods for single-cell RNA-seq data analysis. The Leiden module has a resolution parameter which allows to determine the scale of the partition cluster and therefore the coarseness of the clustering. A higher resolution parameter leads to more clusters.

```python
sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
```

```python
sc.pl.umap(
    adata,
    color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"],
    legend_loc="on data",
)
```

### 11. Annotation
### 12. Data integration

## Inferring trajectories

### 13. Pseudotemporal ordering

Pseudotimes rank cells relative to each other according to their respective stage in the developmental process. Less mature cells are assigned small, mature cells large values.

Probabilistic frameworks assign transition probabilities to ordered cell-cell pairs. Each transition probability quantifies how likely the reference cell is the ancestor of the other cell.

### 14. RNA velocity

Single-cell measurements are snapshot data and can, thus, not be plotted against time. Instead, classical RNA velocity methods rely on studying the cell-specific tuples (u, s) of unspliced and spliced RNA for each gene. The collection of these tuples form the so-called phase portrait. Assuming constant rates of transcription, splicing, and degradation, the phase portraits exhibits an almond shape. The upper arc corresponds to the induction, the lower arc the repression phase. However, as real-world data is noisy, plotting the unspliced against spliced counts does not recover the expected almond shape. Instead, the data needs to be smoothed first. Classically, this preprocessing step consists in averaging the gene expression of each cell over its neighbors in a cell-cell similarity graph.

The first attempt at estimating RNA velocity assumed gene independence and the underlying kinetics to be governed by the above model. Additionally, it is assumed that (1) kinetics reached their equilibrium, (2) rates are constant, and (3) there is a single, common splicing rate across all genes. In the following, we will refer to this model as the **steady-state model** due to the first assumption. The steady-states itself are found in the upper right corner of the phase portrait (induction phase) and its origin (repression phase). Based on these extreme quantiles, the steady-state model estimates the steady-state ratio with a linear regression fit. RNA velocity is then defined as the residual to this fit.

Even though the steady-state model can successfully recover the developmental direction in some systems, it is inherently limited by its model assumptions. The two assumptions readily violated are the common splicing rate across genes and that the equilibria are observed during the experiment.

To overcome the limitations of the steady-state model, several extensions have been proposed. The so-far most popular one is the **EM model** implemented in scVelo. The EM model no longer assumes that steady-states have been reached or that genes share a common splicing rate. Additionally, all datapoints are used to infer the full set of parameters as well as a gene and cell specific latent time of the splicing model. The algorithm uses an expectation-maximization (EM) framework to estimate parameters. The unobserved variables found in the E-step consist of each cell’s time and state (induction, repression, or steady-state). All other model parameters are inferred during the M-step.

While the EM model no longer relies on key assumptions of the steady-state model and, thus, is more broadly applicable, the inferred RNA velocity may still violate prior biological knowledge. The reason for such failure cases are mainly two-fold: On the one hand, the EM model continues to assume constant rates. Consequently, whenever these assumption does not hold, for example in erythroid maturation, the inference is incorrect.

Stady-state model. 
```python
scv.tl.velocity(adata, mode="deterministic")
```
EM model. In order to calculate RNA velocity with the EM model, the parameters of splicing kinetics need to be inferred first. The inference is taken care of by scVelo’s recover_dynamics function.
```python
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
```

### 15. Lineage tracing

## Dealing with conditions

### 16. Differential gene expression analysis
### 17. Compositional analysis
### 18. Gene set enrichment and pathway analysis
### 19. Perturbation modeling

## Modeling mechanisms

### 20. Gene regulatory networks
### 21. Cell-cell communication

## Deconvolution

### 22. Bulk deconvolution

## Chromatin Accessibility

### 23. Single-cell ATAC sequencing
### 24. Quality Control
### 25. Gene regulatory networks using chromatin accessibility

## Spatial omics

### 26. Single-cell data resolved in space
### 27. Neighborhood analysis
### 28. Spatial domains
### 29. Spatially variable genes
### 30. Spatial deconvolution
### 31. Imputation

## Surface protein

### 32. Quality control
### 33. Normalization
### 34. Doublet detection
### 35. Dimensionality Reduction
### 36. Batch correction
### 37. Annotation

## Adaptive immune receptor repertoire

### 38. Immune Receptor Profiling
### 39. Clonotype analysis
### 43. Specificity analysis
### 44. Integrating AIR and transcriptomics

## Multimodal integration

### 45. Paired integration
### 46. Advanced integration
