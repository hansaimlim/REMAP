#BENCHMARK CODES (in MATLAB)</br>

## Brief usage of each matlab code</br>
* REMAP_auto.m : Automatically optimize parameters and get normalized prediction scores based on the optimal parameters.</br>
* REMAP_csv.m : Calculate True Positive Rate (TPR) at top 1% prediction for a train and test file pair.</br>
* REMAP.m : Get the normalized prediction score matrix (matrix Y) for the given dataset using the default parameters.</br>
* REMAP_opt_p6p7.m : Get the optimal parameters (p6 and p7) based on the given 10-fold cross-validation dataset.</br>
* REMAP_opt_rank_iter.m : Get the optimal parameters (rank and iteration) based on the given 10-fold cross-validation datsaet.</br>
* FindTrues.m and TPRbyRowRank.m are used to calculate True Positive Rate (TPR).</br>

## To optimize parameters and get normalized prediction scores, run Y_norm=REMAP_auto(chem_prot_mat,chem_chem_mat,prot_prot_mat)</br>
* The three inputs are chemical-protein association matrix (chem_prot_mat), chemical-chemical similarity matrix (chem_chem_mat), and protein-protein similarity matrix (prot_prot_mat).</br>
* No need to prepare cross-validation datasets. chem_prot_mat should contain ALL known associations.</br>
* It will do grid search to find optimal rank, p6, p7 based on AUC values from 10-fold cross validation.</br>
* Using the optimized parameters, prediction scores are calculated and normalized.</br>
* Please note that the three inputs are in matlab matrix format (.mat files).</br>

## To test on our benchmark datasets, run REMAP_csv(train_csv, test_csv)</br>
* The two inputs are file paths as string. The .csv files contain one chemical-protein pair per line "chemical_index, protein_index".</br>
* REMAP_csv will load the matrices from .csv file, once paths are given</br>
For example:</br>
```
% Example code below
% matlab interactive window. under REMAP/matlab_codes/ directory
>> REMAP_csv('/REMAP/benchmark/NTNL/train_N2L6to10.csv', '/REMAP/benchmark/NTNL/test_N2L6to10.csv')
```
* This will print out the TPR by top 1% prediction (35th rank) for the given test and training files.</br>
* Please note that this code (REMAP_csv) will load and use our benchmark similarity matrices.

--------

## To get the matrix Y (raw prediction score matrix), </br>
 * prepare R: known chemical-protein associations (1 if known, 0 otherwise)</br>
 * chem_chem_mat: chemical-chemical similarity matrix as explained in the reference paper</br>
 * prot_prot_mat: protein-protein similarity matrix as explained in the reference paper</br>
 * and run Y=REMAP(R, chem_chem_mat, prot_prot_mat)</br>
 * Y is a matrix of normalized prediction score</br>
For example:</br>
```
% Example code below
% matlab interactive window. under REMAP/matlab_codes/ directory
>> load /path/to/matrices/R;
>> load /path/to/matrices/chem_chem_mat;
>> load /path/to/matrices/prot_prot_mat;
>> Y=REMAP(matrix_R, chemical_similarity_matrix_name, protein_similarity_matrix_name);
```
 * Matrix Y has the same dimension as matrix_R. Each row represents a chemical, and each column represents a protein</br>
 * For instance, >>Y(1,10) will print the normalized prediction score for the association between the 1st chemical and the 10th protein.
 * Please note that you need to load matrices first, using >>load /path/to/matrix/matrix command</br>

--------

## REMAP parameter optimization.</br>
 * REMAP_opt_p6p7.m and REMAP_opt_rank_iter.m are written for parameter optimization</br>
 * REMAP_opt_p6p7.m optimizes p_chem (p6) and p_prot (p7)</br>
 * REMAP_opt_rank_iter.m optimizes rank (r) and p_iter (iteration)</br>
 * 10-fold cross validation is used for parameter optimization</br>
 * Prepare 10-folded data, chemical-chemical and protein-protein similarity matrices</br>
 * The input directory must contain 10-folded data, the same format as in '/REMAP/benchmark/cv10/'</br>
 * The train and test file names must meet the format, "trainX_10cv.csv" and "testX_10cv.csv", where X is from 1 to 10.</br>
 * The optimization codes will print the optimal parameters based on the 10-fold cross validation.</br>
```
% Example code below
% matlab interactive window. under REMAP/matlab_codes/ directory
>> load ../benchmark/sim/chem_chem_zinc;
>> load ../benchmark/sim/prot_prot_zinc;

% optimize rank and max_iteration
>> REMAP_opt_rank_iter('../benchmark/cv10/',chem_chem_zinc,prot_prot_zinc)

% optimize p_chem (p6) and p_prot (p7)
>> REMAP_opt_p6p7('../benchmark/cv10/',chem_chem_zinc,prot_prot_zinc)
```
--------

## Functions for calculating True Positive Rate.</br>
 * FindTrues.m and TPRbyRowRank.m are required to calculate True Positive Rate.</br>
 * WeightNormalize.m normalizes the raw prediction scores.</br>
 * These functions are called within other REMAP codes above.</br>

--------


```
REFERENCES
[1] Lim, H., Poleksic, A., Yao, Y., Tong, H., He, D., Zhuang, L., Meng, P. and Xie, L., 2016. Large-Scale Off-Target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing. PLOS Comput Biol, 12(10), p.e1005135.
[2] Yao, Y., Tong, H., Yan, G., Xu, F., Zhang, X., Szymanski, B.K. and Lu, J., 2014, November. Dual-regularized one-class collaborative filtering. In Proceedings of the 23rd ACM International Conference on Conference on Information and Knowledge Management (pp. 759-768). ACM.
```
