
# Tri-factorization version of REMAP

 * Use TriFacREMAP_optimization.m to obtain optimal parameters </br>
 ```
 optParams=TriFacREMAP(chem_prot_observed, chem_chem_similarity, prot_prot_similarity, './optimization_result.txt');
 ```
 * Run TriFacREMAP.m to obtain low-rank matrices, and use low-rank matrices for prediction </br>
 ```
 [U,S,V]=TriFacREMAP(chem_prot_observed, chem_chem_similarity, prot_prot_similarity, optParams);
 Y_hat=U*S*V';     %note that Y_hat includes prediction scores for known pairs as well.
 Y_hat=Y_hat.*(~chem_prot_observed);     %removes prediction scores for known pairs.
 ```

------

```
REFERENCES
[1] Hansaim Lim and Lei Xie. “Target gene prediction of transcription factor using a new neighborhood-regularized tri-factorization one-class collaborative filtering algorithm” Proceedings in the 9th ACM Conference on Bioinformatics, Computational Biology, and Health Informatics (ACM BCB) (2018).
[2] Annie Wang, Hansaim Lim (co-first author), Shu-Yuan Cheng, and Lei Xie. “ANTENNA, a multi-rank, multi-layered recommender system for inferring reliable drug-gene-disease associations: repurposing diazoxide as a targeted anti-cancer therapy.” IEEE/ACM Transactions on Computational Biology and Bioinformatics (2018). DOI: 10.1109/TCBB.2018.2812189
```
