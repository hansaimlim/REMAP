#BENCHMARK CODES WRITTEN in MATLAB</br>

## To test on our benchmark datasets, run REMAP_csv(train_csv, test_csv)</br>
For example:</br>
```
>> REMAP_csv('/REMAP/benchmark/NTNL/train_N2L6to10.csv', '/REMAP/benchmark/NTNL/test_N2L6to10.csv')
```
 * The two inputs are file paths as string. REMAP_csv will load the matrices once paths are given.
--------

## To get the matrix Y (raw prediction scores), </br>
 * prepare R: known chemical-protein associations (1 if known, 0 otherwise)</br>
 * chem_chem_mat: chemical-chemical similarity matrix as explained in the reference paper</br>
 * prot_prot_mat: protein-protein similarity matrix as explained in the reference paper</br>
 * and run Y=REMAP(R, chem_chem_mat, prot_prot_mat)</br>
For example:</br>
```
>> load /path/to/matrices/R;
>> load /path/to/matrices/chem_chem_mat;
>> load /path/to/matrices/prot_prot_mat;
>> Y=REMAP(matrix_R, chemical_similarity_matrix, protein_similarity_matrix);
```
 * Matrix Y has the same dimension as matrix_R. Each row represents a chemical, and each column represents a protein.</br>
 * Please note that you need to load matrices first, using >>load /path/to/matrix/matrix command.

--------

## REMAP parameter optimization.</br>
 * REMAP_opt_p6p7.m and REMAP_opt_rank_iter.m are written for parameter optimization.</br>
 * REMAP_opt_p6p7.m optimizes p_chem (p6) and p_prot (p7).</br>
 * REMAP_opt_rank_iter.m optimizes rank (r) and p_iter (iteration).</br>
 * 10-fold cross validation is used for parameter optimization.
 * Prepare 10-folded data, chemical-chemical and protein-protein similarity matrices.</br>
 * The input directory must contain 10-folded data, the same format as in '/REMAP/benchmark/cv10/'.</br>

--------



```
REFERENCES
[1] Hansaim Lim, Aleksandar Poleksic, Hanghang Tong, Yuan Yao, Di He, Luke Zhuang, Patrick Meng, and Lei Xie, "Large-scale Off-target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing" , under review
[2] Yuan Yao, Hanghang Tong, Guo Yan, Feng Xu, Xiang Zhang, Boleslaw K. Szymanski, and Jian Lu. "Dual-regularized one-class collaborative filtering." In Proceedings of the 23rd ACM International Conference on Conference on Information and Knowledge Management, pp. 759-768. ACM, 2014.
```
