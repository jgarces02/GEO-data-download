# GEO-data-download
Code to download microarray data from the GEO database given the GSE and GPL ids.

Use `Download_GEO_Microarray_Data.py` to download expression levels of a set of genes for all all samples:
```
python Download_GEO_Microarray_Data.py GSEid GPLid Gene_list_file
```
The `Gene_list_file` includes the names of genes of interest, one on each line.

Use the other two python codes to download expression data required for the calculation of an EMT scoring metric as described previously by [George *et al*.](https://cancerres.aacrjournals.org/content/77/22/6415.short):
```
python Download_GEO_Microarray_Data_For_EMT_Score.py GSEid GPLid
```

This code was originally written for use in the analysis shown in [Tripathi *et al*.](https://www.frontiersin.org/articles/10.3389/fonc.2018.00244/full) and later adapted for the calculation of the EMT scoring metric.

When using this code for purposes other than calculation of the EMT scoring metric, please cite:
Tripathi S, Jolly MK, Woodward WA, Levine H, and Deem MW (2018) Analysis of Hierarchical Organization in Gene Expression Networks Reveals Underlying Principles of Collective Tumor Cell Dissemination and Metastatic Aggressiveness of Inflammatory Breast Cancer. Front. Oncol. 8:244. doi: 10.3389/fonc.2018.00244.
When using this code for calculation of the EMT scoring metric, also cite:
George, J. T., Jolly, M. K., Xu, S., Somarelli, J. A., & Levine, H. (2017). Survival Outcomes in Cancer Patients Predicted by a Partial EMT Gene Expression Scoring Metric. Cancer Res., 77(22), 6415-6428. doi: 10.1158/0008-5472.CAN-16-3521