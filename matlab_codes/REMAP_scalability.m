function REMAP_scalability()
k=5; %repeats
maxIte=100;

ranks=[100,200,300,400,500];
load('../data/ZINC_ChEMBL_DrugBank/chem_chem_ZCD.mat');
load('../data/ZINC_ChEMBL_DrugBank/chem_prot_ZCD.mat');
load('../data/ZINC_ChEMBL_DrugBank/prot_prot_ZCD.mat');
times=zeros(numel(ranks),k);
for i=1:numel(ranks)
    rank=ranks(i);
    for j=1:k
        times(i,j)=Time_REMAP(chem_prot_ZCD,chem_chem_ZCD,prot_prot_ZCD,rank,maxIte);
    end
end
times_output=zeros(numel(ranks),2);
times_output(:,1)=mean(times,2);
times_output(:,2)=std(times,1,2); %standard deviation normalized by num_obs
writematrix(times_output,'matREMAP_time_on_ZCD_saturn.csv');
clear chem_prot_ZCD chem_chem_ZCD prot_prot_ZCD;

ranks=[100,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000];
load('../BiDD/matrix/chem_chem_zinc.mat');
load('../BiDD/matrix/prot_prot_zinc.mat');
load('../BiDD/matrix/chem_prot_zinc.mat');
times=zeros(numel(ranks),k);
for i=1:numel(ranks)
    rank=ranks(i);
    for j=1:k
        times(i,j)=Time_REMAP(chem_prot_zinc,chem_chem_zinc,prot_prot_zinc,rank,maxIte);
    end
end
times_output=zeros(numel(ranks),2);
times_output(:,1)=mean(times,2);
times_output(:,2)=std(times,1,2); %standard deviation normalized by num_obs
writematrix(times_output,'matREMAP_time_on_ZINC_saturn.csv');

clear chem_prot_zinc chem_chem_zinc prot_prot_zinc i j ranks;
ranks=[100,200,300,400,500];
load('../data/synthetic/chem_chem_syn.mat');
load('../data/synthetic/chem_prot_syn.mat');
load('../data/synthetic/prot_prot_syn.mat');
times=zeros(numel(ranks),k);
for i=1:numel(ranks)
    rank=ranks(i);
    for j=1:k
        times(i,j)=Time_REMAP(chem_prot_syn,chem_chem_syn,prot_prot_syn,rank,maxIte);
    end
end
times_output=zeros(numel(ranks),2);
times_output(:,1)=mean(times,2);
times_output(:,2)=std(times,1,2); %standard deviation normalized by num_obs
writematrix(times_output,'matREMAP_time_on_synthetic_saturn.csv');

clear chem_prot_zinc chem_chem_zinc prot_prot_zinc i j ranks;

end

function elapsed=Time_REMAP(R, chem_chem_mat, prot_prot_mat,rank, maxIte)
%[1] Hansaim Lim, Aleksandar Poleksic, Hanghang Tong, Yuan Yao, Di He, Luke Zhuang, Patrick Meng, and Lei Xie, "Large-scale Off-target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing" , under review
%chem_chem_mat and prot_prot_mat must be square matrices
%number of rows and number of columns of R must match the number of chemicals and proteins, respectively
%prediction matrix Y does NOT show scores for known associations
%in other words, scores for known association pairs are set to 0
%Please refer to the paper for detail

% R : Known association matrix (m by n)  *Required
% chem_chem_mat : row-side similarity score matrix (m by m) (e.g. chemical-chemical similarity)  *Required
% prot_prot_mat : column-side similarity score matrix (n by n) (e.g. protein-protein similarity) *Required
% outfile : output file  *Required

% Normalization : Use score normalization according to the reference paper. Default option is FALSE.
%                  
% cutoff : cutoff score. (between 0 and 1; default 0.5) Raw score or Normalized score depending on the normalization option.
% para : user-defined parameters. (e.g. [0.1, 0.1, 0.01, 200, 300, 0.75, 0.1])


%get number of chemical and protein
m=size(chem_chem_mat, 1);
n=size(prot_prot_mat, 1);

summ = sum(chem_chem_mat,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_mat;

sumn = sum(prot_prot_mat,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_mat;

clear Dm Dn m n
tic;
[U, V] = updateUV(R, Lu, Lv, rank, maxIte);
elapsed=toc;


end

function [U, V] = updateUV(R, Lu, Lv, rank, maxIte)
[m, n] = size(R);
alpha = 0.1;
w = 0.1;
r = 0.1;

gamma = 0.1;
lambda = 0.1;
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
