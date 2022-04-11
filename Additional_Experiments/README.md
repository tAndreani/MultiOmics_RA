### Flow cytometry analysis with R
Raw data and analysis results are stored in: s3://snf-mgln-immunometabolism/CIA_FH2000/FACS

The scripts are:
https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Snippet_Step0_QC.R
https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Snippet_Step1_Read_Comp_Trans.R
https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Snippet_Step2_autoGating_CD45cells.R
https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Snippet_Step3_CATALYST.R
https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Snippet_Step4_ridgelinePlot.R
https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Snippet_Step5_diffExpression.R
https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Snippet_manualGating.R

### Seahorse Real-Time Cell Metabolic Analysis
R Markdown interface on Seahorse raw data: https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Seahorse_11242020.Rmd

The following figures are produced by script: https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Seahorse_Calculation_20210110.R

![Seahorse experiment scheme](https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Demo_OCR_ECAR_20210110.png)

![Mean and SD of raw data on each Seahorse measurement](https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Mean_sd_20210110.png)

![Violin plot on each Seahorse measurement](https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Calculation_20210110.png)

![Dot plot ECAR vs. OCR: Basal condition](https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Basal_OCR_ECAR_20210110.png)

![Dot plot ECAR vs. OCR: Maximal condition](https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Maximal_OCR_ECAR_20210110.png)

The following figures are produced by script: https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Seahorse_Violin_Plot_20210112.R

![Violin plot OCR/ECAR ratio](https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/OCR_ECAR_Ratio_Label_With_Score_20210118.png)

![Violin plot: Basal respiration and glycolysis](https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Calculation_Label_With_Score_Dept_Meeting_20210903.png)

![Violin plot for TA review](https://github.com/tAndreani/MultiOmics_RA/blob/main/Additional_Experiments/Calculation_Label_With_Score_TA_Review_2panels_20210511.png)
