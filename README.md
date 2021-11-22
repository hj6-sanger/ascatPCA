# ascatPCA

ascatPCA requires a panel of normal dataset for denoising copy number variation profiles using Princinpal Component Analysis.


cal_log.py : Calculation of logR based on allele.count files
 Performing PCA with case logR and control logR matrix
*  Input : metadata (e.g., pca_matadata.txt) & logR values from STEP1
*  output : denoised tumor logR data for ASCAT
