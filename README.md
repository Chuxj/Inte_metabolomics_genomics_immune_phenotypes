# Inte_metabolomics_genomics_immune_phenotypes

Codes for generating the major results reported in the manuscript "Integration of metabolomics, genomics and immune phenotypes reveals the causal roles of metabolites in disease" by Chu et al.

1. Cor_analysis: An R code for correlation analsyis between metabolites and immune phenotypes

2. cortest_FADS2_exp_cytokine.R: An R code for correlation analysis between FADS2 gene expression and cytokine level

3. qtl_mapping.R: An R code for qtl mapping with R package MatrixEQTL

4. MR.R: An R code for Mendelian Randomization

5. make_box_plot.pl and boxplot_function_in_perl.R: Codes for batch processing boxplots.

6. PostProcessing_v2.R, ExtractLocus.R: R codes for postprocessing output file from MatrixEQTL.
It takes txt files generated by R package MatrixEQTL as input and generated locus files which could be further used as input on "http://locuszoom.org" to generate locus zoom plot and be used by Colocolization.r to make coloc analysis.

7. Colocolization.r, ColocProcessing2.r: R codes for colocolization analysis and visualization.

8. prediction_function.R: An R function for predicting cytokine production by elastic net

9. manha.R: An R code for making manhattan plot

