# Multi-omics profiling of collagen-induced arthritis mouse model reveals early metabolic dysregulation via SIRT1 axis

Rheumatoid arthritis (RA) is characterized by joint infiltration of immune cells and synovial inflammation which leads to progressive disability. Current treatments improve the disease outcome, but the unmet medical need is still high. New discoveries over the last decade have revealed the major impact of cellular metabolism on immune cell functions. So far, a comprehensive understanding of metabolic changes during disease development, especially in the diseased microenvironment, is still limited. Therefore, we studied the longitudinal metabolic changes during the development of murine arthritis integrating metabolomics and bulk RNA-seq data. We identified an early change in macrophage pathways which was accompanied by oxidative stress, a drop in NAD+ level and induction of glucose transporters. We discovered inhibition of SIRT1, a NAD-dependent histone deacetylase and confirmed its dysregulation in human macrophages and synovial tissue of RA patients. Mining this database should enable the discovery of novel metabolic targets and therapy opportunities in RA. 


##### Check our preprint [Multi-omics profiling of collagen-induced arthritis mouse model reveals early metabolic dysregulation via SIRT1 axis](https://www.biorxiv.org/content/10.1101/2022.03.09.483621v1) in _bioRxiv_.

### Figure 1: Arthritis progression in CIA mice.
#### [Fig.1c Evolution of arthritis severity by scores (sum of 4 paws, plotted as mean ± SEM).](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/linePlot_scoreSum.R)

#### [Fig.1d Arthritis severity at each sacrifice time point per mouse subject.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/boxPlot_scoreSum.R)

#### [Fig.1e Spleen weight measured at each sacrifice time point.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/boxPlot_Spleen_Weight.R)

#### [Fig.1f Plasma TNF-α concentration.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/boxPlot_Plasma_TNF.R)

#### [Fig.1g Plasma CXCL1 concentration.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/boxPlot_Plasma_CXCL1.R)

#### [Fig.1h Plasma CCL5 concentration.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/boxPlot_Plasma_CCL5.R)

### Figure 2: UPLC-MS/MS metabolomics revealed early metabolic alterations in plasma of CIA mice.

#### [Fig.2a Principal component analysis (PCA) score plot with plasma samples projected onto the first 2 principal components.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/PCA_Plasma_Metabolome_Manuscript.R)

#### [Fig.2b Percentage of each metabolic group detected by UPLC-MS/MS.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/donutChart_Manuscript.R)

#### [Fig.2c Number of significant metabolites at each time point (CIA samples compared to the controls).](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/barPlot_DEG_Met_Number.R)

#### [Fig.2d Percentage of significant metabolites per metabolite group at each time point.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/linePlot_barPlot_metabolite_proportion.R)

#### [Fig.2e Percentage of upregulated and downregulated significant metabolites per group over time.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/linePlot_barPlot_metabolite_proportion.R)

#### [Fig.2f Fold changes and significances of peptides that are significant at >= 2 time points.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/dotPlot_peptide.R)

#### [Fig.2g Fold changes and significances of energy metabolites and related cofactors that are significant at >= 2 time points.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/dotPlot_carbohydrate.R)

### Figure 3: Transcriptional profiling of CIA paws over disease progression revealed innate immune cell activation and metabolic adaptations.

#### [Fig.3a PCA score plot with paw samples color-coded by their respective arthritis scores.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/PCA_DESeq2_manuscript.R)

#### [Fig.3b PCA score plot with paw samples color-coded by their respective sacrifice time points.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/PCA_DESeq2_manuscript.R)

#### [Fig.3c Number of differentially expressed genes at each time point (CIA samples compared to the controls).](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/barPlot_DEG_Number.R)

#### [Fig.3d Enriched or depleted Gene Ontology biological process gene sets at each time point.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/dotPlot_GSEA_go_bp.R)

#### [Fig.3e Top 10 activated and top 10 inhibited upstream regulators from analysis by IPA.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/heatmap_upstream_regulator_manuscript.R)

### Figure 4: Tissue MALDI-MS imaging showed reduction of NAD+ in CIA paws.

#### [Fig.4a PCA score plot with paw samples color-coded by their respective arthritis scores.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/sum_mz_PCA_mSet.R)

#### [Fig.4b PCA score plot with paw samples color-coded by their respective sacrifice time points.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/sum_mz_PCA_mSet.R)

#### [Fig.4d Image signal quantification of NAD+ on paws.](https://github.com/tAndreani/MultiOmics_RA/blob/main/Codes/boxPlot_NAD_manuscript.R)

### Figure 5: Multi-omics data integration identifies ROS as disease-correlated factor.

#### Fig.5b Proportion of total variance (Var.) explained by individual factors for each omics experiment. 

#### Fig.5c Visualization of samples using Factors 1, 2 and 3. 

#### Fig.5d Absolute weights of the top features of Factors 1 in transcriptome data. 

#### Fig.5e Gene set enrichment analysis on the feature weights of mRNA in Factor 1 (FDR < 0.05).

### Figure 6: Transcriptional profiling of LPS-stimulated human macrophages partially overlaps to early transcript changes in mouse CIA.

#### Fig.6a PCA score plot with four samples of human macrophages before and after LPS stimulation. 

#### Fig.6b Overlap of significant differentially expressed genes between CIA paws at week three and LPS-stimulated macrophages. 

#### Fig.6c Enriched or depleted Gene Ontology biological process gene sets. 

#### Fig.6d Top 20 activated and top 20 inhibited upstream regulators from analysis by IPA. 

### Figure 7: Gene expression of human synovium over disease progression to RA showed similar changes in oxidative stress and innate immunity.

#### Fig.7b PCA score plot with patient samples color-coded by disease stage. 

#### Fig.7c Number of significant differentially expressed genes at each disease stage. 

#### Fig.7d Overlap of significant differentially upregulated genes between CIA paws at week three and human synovium of different RA stages. 

#### Fig.7e Enriched or depleted Gene Ontology biological process gene sets. 

#### Fig.7f Top 10 activated and top 10 inhibited upstream regulators from analysis by IPA.

### Figure 8: Metabolic reprogramming in early RA.

#### Fig.8a-c Differential gene expression of major metabolic genes in the RNA-seq studies of CIA mouse paws, LPS-stimulated human macrophages, and human RA synovium biopsies. 
