function REMAP(R, chem_chem_mat, prot_prot_mat,outfile, Normalization,cutoff,para)
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



maxNumCompThreads(3); %determine the maximum number of cores to use
if nargin < 4
    msg='Too few arguments are given. At least four arguments required.';
    error(msg)
elseif nargin < 5
    %four arguments given
    para = [0.1, 0.1, 0.01, 200, 300, 0.75, 0.1]; %default parameters
    cutoff=0.5;
    Normalization='no';
elseif nargin < 6
    %five arguments given
    para = [0.1, 0.1, 0.01, 200, 300, 0.75, 0.1]; %default parameters
    cutoff=0.5;
elseif nargin < 7
    para = [0.1, 0.1, 0.01, 200, 300, 0.75, 0.1];	% para: p_reg, squared p_weight, p_imp, rank, p_iter, p_chem, p_prot
end
if length(para) < 7
    msg=['The parameter must contain 7 values.\n' ...
        'e.g.)para=[0.1, 0.1, 0.01, 200, 300, 0.75, 0.1]'];
    error(msg)
end
if strcmpi(Normalization,'true')||strcmpi(Normalization,'normalization')||strcmpi(Normalization,'yes')
    msg=['Normalization True, cutoff Normalized score ' num2str(cutoff)];
    disp(msg);
    Normalize=1
elseif strcmpi(Normalization,'false')||strcmpi(Normalization,'raw')||strcmpi(Normalization,'no')
    msg=['Normalization False, cutoff Raw score ' num2str(cutoff)];
    disp(msg);
    Normalize=0
else
    msg=['Normalization False, cutoff Raw score ' num2str(cutoff)];
    disp(msg);
    Normalize=0
end
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
[U, V] = updateUV(R, Lu, Lv, para);
if Normalize
    Y=WeightNormalize(U*V',R,cutoff);
else
    Y=U*V';
end
clear U V;

Y=Y-Y.*R; %remove known associations
[i,j,k]=find(Y);
idxs=horzcat(i,j,k);
dlmwrite(outfile,idxs,'delimiter','\t','precision',10);
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
