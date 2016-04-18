function REMAP_csv(train_csv, test_csv)
%[1] Hansaim Lim, Aleksandar Poleksic, Hanghang Tong, Yuan Yao, Di He, Luke Zhuang, Patrick Meng, and Lei Xie, "Large-scale Off-target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing" , under review
tic;
para = [0.1, 0.1, 0.01, 200, 300, 0.75, 0.1];	% para: p_reg, squared p_weight, p_imp, rank, p_iter, p_chem, p_prot

%chem_chem_zinc and protein_protein_zinc_blast matrices from chem-chem and prot-prot files
load ../benchmark/sim/chem_chem_zinc;
load ../benchmark/sim/protein_protein_zinc_blast;
%get number of chemical and protein
m=size(chem_chem_zinc, 1);
n=size(protein_protein_zinc_blast, 1);
%convert csv to matrix
train_line = csvread(train_csv);
train = sparse(train_line(:,1), train_line(:,2), 1, m, n);      %12384 chemicals and 3500 proteins in ZINC
test = csvread(test_csv);
testmat = sparse(test(:,1), test(:,2), 1, m, n);	%matrix form for test set

summ = sum(chem_chem_zinc,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_zinc;

sumn = sum(protein_protein_zinc_blast,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - protein_protein_zinc_blast;

[U, V] = updateUV(train, Lu, Lv, para);


test_result = TPRbyRowRank(FindTrues(U*V', testmat), 100);   %max cutoff rank 100
disp(['Rank=' num2str(rank) 'Iter=' num2str(Ite) 'TPR at top 35th=' num2str(test_result(35,2)) ])
clear train;
clear test;
clear test_result;
toc

end

function [U, V] = updateUV(R, Lu, Lv, para)
[m, n] = size(R);
alpha = para(1);
w = para(2);
r = para(3);
rank = para(4);
maxIte = para(5);
gamma = para(6);
lambda = para(7);
ite = 0;

U0 = rand(m, rank);
V0 = rand(n, rank);

Lu_plus = (abs(Lu) + Lu) / 2;
Lu_minus = (abs(Lu) - Lu) / 2;

Lv_plus = (abs(Lv) + Lv) / 2;
Lv_minus = (abs(Lv) - Lv) / 2;

while ite <maxIte 

    UVT = get_UVT(R, U0, V0);
    U0 = updateU(R, UVT, w, r, Lu_plus, Lu_minus, U0, V0, alpha, gamma);
    V0 = updateU(R', UVT',w, r, Lv_plus, Lv_minus, V0, U0, alpha, lambda);
    
    ite = ite + 1;
end

U = U0;
V = V0;

end

function [UVT] = get_UVT(R, U, V)

[m, n] = size(R);

[I, J] = find(R);
iSize = size(I, 1);
K = ones(iSize,1);
count = 0;
for j = 1:n
    id = find(R(:,j));
    len = length(id);
    K(count+1:count+len) = U(id, :) * V(j,:)';
    count = count + len;
end

UVT = sparse(I, J, K, m, n);

end

function [U1] = updateU(R, UVT, w, p, Lu_plus, Lu_minus, U0, V, lambda, gamma)

[m,n]=size(R);
U1 = U0 .* sqrt( ((1-w*p)*R*V + ones(m,1)*p*((w*ones(1,n))*V) + gamma .* Lu_minus * U0) ./ ((1-w)*UVT*V + w * (U0*(V'*V))  + gamma.* Lu_plus * U0 + lambda * U0) );

U1(isnan(U1)) = 0;

end
