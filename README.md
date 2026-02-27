# MS_day84
High-Dose Vitamin D3 Supplementation Associates with Innate Immune Chromatin Remodeling and Context-Dependent Transcriptomic Responses in Multiple Sclerosis


The raw data can be found on Geo with accsession no. GSE320190 (RNA seq MS) and GSE320191 (ATAC seq MS), GSE319989 (RNA seq Healthy day 84), GSE278885 (RNA seq Healthy day 0) and GSE303320 (ATAC-seq Healthy day 0).
The supplementary table TABLE S2.xlsx and Table S3.xlsx (link on the paper) have normalised counts from rna seq data and atac seq data respectively of both MS and healthy individuals.

Vitamin D deficiency is a recognized risk factor for multiple sclerosis (MS), yet the in vivo molecular mechanisms linking supplementation to immune regulation remain poorly defined. Here, we performed longitudinal RNA-seq and ATAC-seq profiling of peripheral blood mononuclar cells (PBMCs) from MS patients who either initiated daily high-dose vitamin D3 dosing or continued regular supplementation, alongside healthy individuals receiving monthly bolus supplementation. Across three months, MS samples remained transcriptionally and epigenetically distinct from healthy controls. Although cohort heterogeneity precludes causal inference, the longitudinal design enabled within-individual assessment of regimen-associated molecular changes. Vitamin D3 supplementation did not induce global transcriptome-wide reprogramming but instead modulated discrete, context-dependent immune gene networks. In MS patients receiving daily high-dose dosing, vitamin D-responsive genes were enriched for innate immune and inflammatory pathways, including Toll-like receptor and cytokine-regulatory processes. Chromatin accessibility profiling revealed widespread MS-associated regulatory differences and partial overlap between regulatory and transcriptional signals. Intersecting vitamin D-responsive genes with published vitamin D receptor binding data highlighted a subset of candidates with plausible MS relevance, motivating targeted mechanistic follow-up. Together, these findings define regimen-specific, context-dependent immune regulatory programs associated with vitamin D₃ supplementation in MS and demonstrate the added value of integrating chromatin accessibility profiling in real-world nutritional interventions


The file 25OHD_levels.R has code for 25OHD level comparision for all 3 groups of high dose MS , regular dose MS and Healthy individuals taking vitamin D bolus for day 0 and day 84.

The file ATAC_vdr_rna_accesibility and distance from tss.R shows the integration of RNA seq and ATAC seq along with previously published VDR chip seq data

The file MA_DGE.R has the code for generating MA plots and differential gene expresion for all the 3 groups.

The file pheatmap_MS_Healthy.R has the code for pheatmap (heirarchical clustering) showing the overall transcriptomic difference between MS patients and healthy individuals pheatmap.
The same code has been used to make pheatmap for differentially expressed genes betwwen MS and healthy for day 0 i.e baseline expression differences using 1480 genes which can be filtered from Table S2.
The similar approach has been used for generating pheatmap for ATAC seq data (k means clustering)

