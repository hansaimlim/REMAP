
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
