# Cubic-Spline-based-Temporal-Analysis-Workflow

Our cubic spline-based computational workflow is widely applicable. It can capture the distinct temporal molecular signatures associated with the temporal profile of various types of disease pathogenesis. 

Cubic spline-based temporal clustering:
Mal-biotin labeled cysteine sites or total cysteine sites with abundance value for at least 4 out of 6 time points in both ISO and Vehicle groups were selected for cubic spline-based temporal clustering. The averaged ratio of abundance in ISO to Vehicle across replicates was calculated per site. Missing abundance values for a modification site were imputed using average abundance of remaining time-points. After scaling and centering, cubic splines were fitted to ratios across the time points in the R statistical programming language (v3.4.3). The predicted abundance ratios from cubic spline were used for K-mean clustering (kmean package in R) following identification of the cluster numbers (mclust package in R). 
