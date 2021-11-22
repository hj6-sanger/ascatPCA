# ascatPCA

# ascatPca

ASCAT for metanorm style analysis requiring Panel of Normal data for PCA processing.


STEP 1 : Calculation of logR based on allele.count file (cal_log.py). 

*  Input : list of matched normal and tumor (e.g., ascat_sample.txt) & count file from ASCAT allele counter
*  Output :logR values for each tumor.

STEP 2 : Performing PCA with case logR and control logR matrix
*  Input : metadata (e.g., pca_matadata.txt) & logR values from STEP1
*  output : denoised tumor logR data for ASCAT
