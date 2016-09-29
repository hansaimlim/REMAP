function REMAP_opt_rank_p6_p7(input_dir,chem_chem_mat,prot_prot_mat, iter)
%[1] Hansaim Lim, Aleksandar Poleksic, Hanghang Tong, Yuan Yao, Di He, Luke Zhuang, Patrick Meng, and Lei Xie, "Large-scale Off-target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing" , under review

%10fold cross validation for ZINC dataset
%input_dir='../benchmark/cv10/';
%chem_chem_zinc and protein_protein_zinc_blast matrices from chem-chem and prot-prot files
%load ../benchmark/sim/chem_chem_zinc;
%load ../benchmark/sim/protein_protein_zinc_blast;
output_file='./REMAP_rank_p6_p7_optimization.txt';
ranks = [100, 200, 300, 400, 500]; %rank parameters to be optimized. Should NOT exceed the smallest matrix dimension
cutoff_rank=35;
%get number of chemical and protein
m=size(chem_chem_mat, 1);
n=size(prot_prot_mat, 1);

if min(m,n) < max(ranks)
   msg=['Rank parameter cannot exceed the smallest matrix dimension(' num2str(min(m,n)) '). Please correct the rank parameter range to optimize.'];
   error(msg);
end
summ = sum(chem_chem_mat,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_mat;

sumn = sum(prot_prot_mat,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_mat;

disp(['Parameter optimization Start'])
disp(['rank, p_chem(p6), and p_prot(p7) are being optimized using cross validation'])
disp(['iteration=' num2str(iter)])
disp(['output=' output_file])


p6s = 0:0.1:1;
p7s = 0:0.1:1;

best_paras=zeros(1,4); %[TPR rank p6 p7] at best TPR
tpr_best=0;
fileid=fopen(output_file,'a+t');
for i = 1:numel(ranks)
	for j = 1:numel(p6s)
        for l = 1:numel(p7s)
            para = [0.1, 0.1, 0.01, ranks(i), iter, p6s(j), p7s(l)]; % para: lambda, squared global weight, r, rank, maxIte, p_chem, p_prot
            TPR35_sum=0;            %sum of TPR35 for each cross validation
            for k=1:10
             trainfile=[input_dir 'train' num2str(k) '_10cv.csv'];
             testfile =[input_dir 'test' num2str(k) '_10cv.csv'];
             trline=csvread(trainfile);
             tsline=csvread(testfile);
             TR=sparse(trline(:,1), trline(:,2), 1, m, n);
             TS=sparse(tsline(:,1), tsline(:,2), 1, m, n);
             [U, V] = updateUV(TR, Lu, Lv, para);
             test_result = TPRbyRowRank(FindTrues(U*V', TS), cutoff_rank);   %max cutoff rank 35
             TPR35_sum = TPR35_sum + test_result(cutoff_rank,2);
             clear U V TR TS trline tsline test_result;
            end
            tpr35=(TPR35_sum/10);
            if tpr35 > best_paras(1,1)
               best_paras(1,:)=[tpr35, ranks(i), p6s(j), p7s(l)]; %update best parameters               
            end
            text=['Top' num2str(cutoff_rank) ' TPR: ' num2str(tpr35) ' iter=' num2str(iter) ' rank=' num2str(ranks(i)) ', p6=' num2str(p6s(j)) ', p7=' num2str(p7s(l)) '\n'];
            fprintf(fileid, text, 'char');
        end
	end
end
fprintf(fileid,['best TPR=' num2str(best_paras(1,1)) ' at rank=' num2str(best_paras(1,2)) ' p6=' num2str(best_paras(1,3)) ' p7=' num2str(best_paras(1,4)) ],'char');
clear train;
clear test;
clear test_result;
disp(['Parameter optimization complete\n'])
fclose(fileid);

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
