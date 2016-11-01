function Y=REMAP_auto(R, chem_chem_mat, prot_prot_mat)
%[1] Hansaim Lim, Aleksandar Poleksic, Hanghang Tong, Yuan Yao, Di He, Luke Zhuang, Patrick Meng, and Lei Xie, "Large-scale Off-target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing" , under review
%This script optimizes parameters using 10-fold cross validation.
%Using the optimized parameters, run REMAP and get the normalized
%prediction score matrix, Y. 

%chem_chem_mat and prot_prot_mat must be square matrices
%number of rows and number of columns of R must match the number of chemicals and proteins, respectively
%Please refer to the paper for detail

maxNumCompThreads(3); %determine the maximum number of cores to use
maxIter=100;
disp('Parameter optimization (grid search) started...\n\r');
best_paras=para_opt(R,chem_chem_mat,prot_prot_mat,maxIter); %parameter optimization
best_auc=best_paras(1);
para = [0.1, 0.1, 0.01, best_paras(2), maxIter, best_paras(3), best_paras(4)];	% para: p_reg, squared p_weight, p_imp, rank, p_iter, p_chem, p_prot
disp('Parameter optimization complete!\n\r');
disp(['Best mean AUC: ' num2str(best_auc) ' at rank=' num2str(best_paras(2)) ', iter=' num2str(maxIter) ', p6=' num2str(best_paras(3)) ', p7=' num2str(best_paras(4))]);
%get number of chemical and protein
m=size(chem_chem_mat, 1);
n=size(prot_prot_mat, 1);

summ = sum(chem_chem_mat,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_mat;

sumn = sum(prot_prot_mat,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_mat;

[U, V] = updateUV(R, Lu, Lv, para);
Y=WeightNormalize(U*V',R); %normalize scores based on the reference
end


function best_paras_auc=para_opt(R,chem_chem_mat,prot_prot_mat,iter)
%[1] Hansaim Lim, Aleksandar Poleksic, Hanghang Tong, Yuan Yao, Di He, Luke Zhuang, Patrick Meng, and Lei Xie, "Large-scale Off-target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing" , under review

%get optimized rank, p6, p7 based on 10-fold cross validation
%return optimized parameters and results
%output_file='./REMAP_rank_p6_p7_optimization.txt';
kFold=10; %10-fold cross-validation for optimization
m=size(chem_chem_mat, 1); %num chem
n=size(prot_prot_mat, 1); %num prot
p6s = 0:0.25:1;
p7s = 0:0.25:1;
if min(m,n)>300
    ranks=100:100:300;
else
    maxrank=min(m,n);
    ranks=floor(maxrank/3):floor(maxrank/3):maxrank;
end

summ = sum(chem_chem_mat,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_mat;

sumn = sum(prot_prot_mat,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_mat;

best_paras_auc=zeros(1,4);
best_paras_aupr=zeros(1,4);
auc_best=0;
aupr_best=0;

gridNumTotal=numel(ranks)*numel(p6s)*numel(p7s); %total number of grid
gridNumComplete=0;
for i = 1:numel(ranks)
	for j = 1:numel(p6s)
        for l = 1:numel(p7s)
            para = [0.1, 0.1, 0.01, ranks(i), iter, p6s(j), p7s(l)]; % para: lambda, squared global weight, r, rank, maxIte, p_chem, p_prot
            %prepare CV object
            positive_idx=find(R>0);
            cv=cvpartition(length(positive_idx),'Kfold',kFold);
            AUC=zeros(1,kFold);
            AUPR=zeros(1,kFold);
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
                 [~,~,~,auc]=perfcurve(TestCompressed(:),P(:),1,'Options',statset('UseParallel',true));
                 [~,~,~,aupr]=perfcurve(TestCompressed(:), P(:), 1, 'Options',statset('UseParallel',true), 'xCrit', 'reca', 'yCrit', 'prec');
                AUC(1,k)=auc;
                AUPR(1,k)=aupr;
            end
            if mean(AUC)>auc_best
                auc_best=mean(AUC); %update best AUC
                best_paras_auc=[auc_best,ranks(i),p6s(j),p7s(l)];
            end
            if mean(AUPR)>aupr_best
                aupr_best=mean(AUPR); %update best AUPR
                best_paras_aupr=[aupr_best,ranks(i),p6s(j),p7s(l)];
            end
            gridNumComplete=gridNumComplete+1;
            disp([num2str(gridNumComplete) ' grid searched. (Total ' num2str(gridNumTotal) ' grid)']);
        end
	end
end

clear Train Test TestCompressed AUC AUPR i j k l para P U V testidx positive_idx;

end

function [U, V] = updateUV(R, Lu, Lv, para)
% para: lambda, r, T, rank, maxIte, ite_of_bisection method, topN
[m, n] = size(R);
alpha = para(1);
w = para(2);
r = para(3);
rank = para(4);
maxIte = para(5);
gamma = para(6);
lambda = para(7);
ite = 0;

%IU = ones(m,n) - R;
%W = IU * w + R;
%IU = sparse(IU) * r;

U0 = rand(m, rank);
V0 = rand(n, rank);

Lu_plus = (abs(Lu) + Lu) / 2;
Lu_minus = (abs(Lu) - Lu) / 2;

Lv_plus = (abs(Lv) + Lv) / 2;
Lv_minus = (abs(Lv) - Lv) / 2;

while ite <maxIte 
    %[RMSE, MAE] = get_diff(test, U0', V0');
    %fprintf('Ite = %d, RMSE = %0.4f, MAE = %0.4f, time = %0.4f\n', ite, RMSE, MAE, etime(clock, t0));
    %[MAP, MPR, HLU, AUC, avgF, avgP, avgR] = get_diff(test, U0, V0, para);
    %fprintf('Ite = %d, MAP = %0.4f, MPR = %0.4f, HLU = %0.4f, AUC = %0.4f, avgF  = %0.4f, avgP = %0.4f, avgR = %0.4f, time = %0.4f\n', ite, MAP, MPR, HLU, AUC, avgF, avgP, avgR, etime(clock, t0));
    

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
