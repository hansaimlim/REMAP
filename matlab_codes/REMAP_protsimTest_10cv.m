function [mean_auc_diff,pval]=REMAP_protsimTest_10cv(R,chem_chem,prot_prot,ProtSimType,ProtMatType)
%[1] Hansaim Lim, Aleksandar Poleksic, Hanghang Tong, Yuan Yao, Di He, Luke Zhuang, Patrick Meng, and Lei Xie, "Large-scale Off-target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing" , under review
%This script is to test different protein-protein similarity metrics
%First using the given training matrix, get the low-rank matrix of
%protein-side (V). Then calculate similarity measurement on the V matrix.
%Run prediction again using the new similarity matrix to see improvement.
%Note that the V must be derived using training dataset only.

maxNumCompThreads(3); %determine the maximum number of cores to use
tic;
MatTypeList={'blank','lin','lin_cut05','path','path_cut05','resnik','resnik_cut05','vector','vector_cut05'};
if ismember(ProtMatType,MatTypeList)
    disp(['prot-prot-sim: ' ProtMatType]);
else
    msg='Please check the similarity matrix type (e.g. renik, lin or blank)';
    error(msg);
end
if strcmp(ProtSimType,'euc')
    ProtSimType='euclidean';
elseif strcmp(ProtSimType,'minkow')
    ProtSimType='minkowski';
elseif strcmp(ProtSimType,'cos')
    ProtSimType='cosine';
elseif strcmp(ProtSimType,'ham')
    ProtSimType='hamming';
elseif strcmp(ProtSimType,'jac')
    ProtSimType='jaccard';
elseif strcmp(ProtSimType,'cor')
    ProtSimType='correlation';
elseif strcmp(ProtSimType,'spear')
    ProtSimType='spearman';
end
SimTypeList={'euclidean','minkowski','cosine','hamming','jaccard','correlation','spearman'}; %available types of similarity
if ismember(ProtSimType,SimTypeList)
    disp(['testing prot-prot sim: ' ProtSimType]);
else
    msg='Please check the similarity type.\n';
    error(msg);
end


%calculate AUC based on 10 fold cross validation
%defalut parameters below. adjust low-rank based on the size
para = [0.1, 0.1, 0.01, 200, 300, 0.75, 0.1];	% para: p_reg, squared p_weight, p_imp, rank, p_iter, p_chem, p_prot
kFold=10;

%get number of chemical and protein
m=size(chem_chem, 1);
n=size(prot_prot, 1);
if min(m,n) < 200
    para(4)=floor(min(m,n)/5);
	disp(['Rank parameter adjusted to r=' num2str(para(4))]);
end
%convert csv to matrix

summ = sum(chem_chem,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem;

sumn = sum(prot_prot,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot;


%prepare CV object
positive_idx=find(R>0);
cv=cvpartition(length(positive_idx),'Kfold',kFold);
AUC=zeros(1,kFold);
AUC_simtype=zeros(1,kFold);
%AUPR=zeros(1,kFold);
for k=1:kFold
	%Cross validation, K-fold
	testidx=positive_idx(test(cv,k)); %k-th test set
	Train=R;
	Train(testidx)=0; %mask test pairs to 0
	Test=R-Train;
	TestCompressed=Test(sum(Test,2)>0,:);%only rows with nonzero element
	[U, V] = updateUV(Train, Lu, Lv, para);
	P=U*V';
	P=P(sum(Test,2)>0,:); %same rows as test compressed matrix
%    [auc,aupr]=AUROC(TestCompressed,P);
	 [~,~,~,auc]=perfcurve(TestCompressed(:),P(:),1,'Options',statset('UseParallel',true));%auc for default test
%	 [~,~,~,aupr]=perfcurve(TestCompressed(:), P(:), 1, 'Options',statset('UseParallel',true), 'xCrit', 'reca', 'yCrit', 'prec');
	AUC(1,k)=auc; %update auc for default test
    
    %start protsime test
    prot_prot_new=squareform(pdist(V,ProtSimType)); 
    prot_prot_new=1-prot_prot_new./max(prot_prot_new(:));
    sumn_new=sum(prot_prot_new,2);
    Dn_new=spdiags(sumn_new,0,n,n);
    Lv_new=Dn_new-prot_prot_new;
    [U_new,V_new]=updateUV(Train,Lu,Lv_new,para);
    P_new=U_new*V_new';
    P_new=P_new(sum(Test,2)>0,:);
    [~,~,~,auc_new]=perfcurve(TestCompressed(:),P_new(:),1,'Options',statset('UseParallel',true)); %auc for new prot-sim
    AUC_simtype(1,k)=auc_new;
%	AUPR(1,k)=aupr;
end

mean_AUC=mean(AUC);
std_AUC=std(AUC);
mean_AUC_new=mean(AUC_simtype);
std_AUC_new=std(AUC_simtype);
mean_auc_diff=mean_AUC_new-mean_AUC;
% std_auc_diff=std_AUC_new-std_AUC;
[~,pval,~,~]=ttest2(AUC,AUC_simtype);

% disp(['mean AUC(default)=' num2str(mean_AUC) ', std. AUC(default)=' num2str(std_AUC) ' at r=' num2str(para(4)) ' kFold=' num2str(kFold)]);
% disp(['mean AUC(' ProtSimType ')=' num2str(mean_AUC_new) ', std. AUC(' ProtSimType ')=' num2str(std_AUC_new) ' two-sample t-test p-value=' num2str(pval)]);
% disp(['mean AUC improvement=' num2str(mean_auc_diff) ', p-value=' num2str(pval)]);
% disp(ProtSimType);
output=[mean_AUC, std_AUC, mean_AUC_new, std_AUC_new, mean_auc_diff, pval];
disp(output);
toc
clear P U V AUC Train Test;

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
