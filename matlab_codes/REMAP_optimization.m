function REMAP_optimization(R,chem_chem,prot_prot)
%[1] Lim, H., Poleksic, A., Yao, Y., Tong, H., He, D., Zhuang, L., Meng, P. and Xie, L., 2016. Large-Scale Off-Target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing. PLOS Comput Biol, 12(10), p.e1005135.
maxNumCompThreads(3); %determine the maximum number of cores to use
%calculate AUC based on 10 fold cross validation
%defalut parameters below. adjust low-rank based on the size
kFold=10; %10-fold cross validation for optimization
fid=fopen('../REMAP_optimization_ZINC_quick.txt','at+');
%get number of chemical and protein
m=size(chem_chem, 1);
n=size(prot_prot, 1);
% if min(m,n) < 200
%     para(4)=floor(min(m,n)/5);
% 	disp(['Rank parameter adjusted to r=' num2str(para(4))]);
% end

summ = sum(chem_chem,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem;

sumn = sum(prot_prot,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot;

ranks=100:100:300;
% iters=100:100:200;
iter=100;
p6s=0:0.33:1.0;
p7s=0:0.33:1.0;

%prepare CV object
positive_idx=find(R>0);
cv=cvpartition(length(positive_idx),'Kfold',kFold);
best_paras=zeros(1,4);
best_auc=0;
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n','LowRank','MaxIter','p6','p7','avgAUC','stdAUC');
for ir=1:numel(ranks)
%     for ii=1:numel(iters)
        for i6=1:numel(p6s)
            for i7=1:numel(p7s)
                para = [0.1, 0.1, 0.01, ranks(ir), iter, p6s(i6), p7s(i7)];
                AUC=zeros(1,kFold);
                tic;
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
                    [~,~,~,auc]=perfcurve(TestCompressed(:),P(:),1);
                    AUC(1,k)=auc;
                end
                avgauc=mean(AUC);
                fprintf(fid,'%5d\t%5d\t%.3g\t%.3g\t%.5g\t%.5g\n',ranks(ir),iter,p6s(i6),p7s(i7),avgauc,std(AUC));
                if avgauc>best_auc
                    best_auc=avgauc;
                    best_paras=[ranks(ir),iter,p6s(i6),p7s(i7)];        
                end
                toc
            end
        end
%     end
end

msg=['Best avg. AUC=' num2str(best_auc) ' at rank=' num2str(best_paras(1)) ', iter=' ...
    num2str(best_paras(2)) ', p6=' num2str(best_paras(3)) ', p7=' num2str(best_paras(4)) ];
disp(msg);
fprintf(fid,'%s\n',msg);
fclose(fid);
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
