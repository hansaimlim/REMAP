function TriFacREMAP_PerfEval(R, chem_chem_mat, prot_prot_mat, param, outfile)
%evaluate performance of tREMAP using given parameters
%maxNumCompThreads(2); %determine the maximum number of cores to use
kFold = 10;
fid=fopen(outfile,'at+');
lowrank_u=param(1); %chem-side low rank
lowrank_v=param(2); %prot-side low rank
p_chem=param(3); %sometimes called p6
p_prot=param(4); %sometimes called p7
para=[0.1,0.1,0.01,lowrank_u,lowrank_v,200,p_chem,p_prot];
if length(para) < 8
    msg=['The parameter must contain 8 values.\n' ...
        'e.g.)para=[0.1, 0.1, 0.01, 200, 50, 300, 0.75, 0.1]'];
    error(msg)
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

%prepare CV object
positive_idx=find(R>0);
cv=cvpartition(length(positive_idx),'Kfold',kFold);

    topKrec=zeros(kFold,n);
    AUC=zeros(1,kFold);
    MAP=zeros(1,kFold); HLU=zeros(1,kFold); MPR=zeros(1,kFold);
    for k=1:kFold
        tic;
        %Cross validation, K-fold
        testidx=positive_idx(test(cv,k)); %k-th test set
        Train=R;
        Train(testidx)=0; %mask test pairs to 0
        Test=R-Train;
        TestCompressed=Test(sum(Test,2)>0,:);%only rows with nonzero element
        [U, S, V] = updateUSV(Train, Lu, Lv, para);
        Pred=U*S*V';
        Pred=Pred(sum(Test,2)>0,:); %same rows as test compressed matrix
        toc
        topkrec=TopKRecall(Pred,TestCompressed,350);
        topKrec(k,1:size(topkrec,2))=topkrec(2,:);
        [~,~,~,auc]=perfcurve(TestCompressed(:),Pred(:),1);     
        [map,hlu,mpr]=evaluatePerf(TestCompressed,Pred);
        AUC(1,k)=auc; MAP(1,k)=map; HLU(1,k)=hlu; MPR(1,k)=mpr;
    end
    avgAUC=mean(AUC);stdAUC=std(AUC);avgHLU=mean(HLU);stdHLU=std(HLU);
    avgMAP=mean(MAP);stdMAP=std(MAP);avgMPR=mean(MPR);stdMPR=std(MPR);
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','AvgAUC','StdAUC','AvgHLU','StdHLU',...
        'AvgMAP','StdMAP','AvgMPR','StdMPR');
    fprintf(fid,'%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n',avgAUC,stdAUC,...
        avgHLU,stdHLU,avgMAP,stdMAP,avgMPR,stdMPR);
    fprintf(fid,'%s\t%s\t%s\n','K','AvgTopKRecall','StdTopKRecall');
    avgTopKRec=mean(topKrec);
    stdTopKRec=std(topKrec);
    for i=1:size(avgTopKRec,2)
        fprintf(fid,'%3d\t%.5g\t%.5g\n',i,avgTopKRec(i),stdTopKRec(i));
    end
fclose(fid);
end

function [U, S, V] = updateUSV(R, Lu, Lv, para)
[m, n] = size(R);
alpha = para(1);
w = para(2);
%W = R+w*(~R); %weight as matrix
r = para(3);
%P = r*(~R); %imputation as matrix
rank_u = para(4);
rank_v = para(5);
maxIte = para(6);
gamma = para(7);
lambda = para(8);
ite = 0;

U0 = rand(m, rank_u);
V0 = rand(n, rank_v);
S0 = rand(rank_u, rank_v);
Lu_plus = (abs(Lu) + Lu) / 2;
Lu_minus = (abs(Lu) - Lu) / 2;

Lv_plus = (abs(Lv) + Lv) / 2;
Lv_minus = (abs(Lv) - Lv) / 2;

while ite <maxIte 

    UVT = get_UVT(R, U0, S0, V0);
    U0 = updateU(R, UVT, w, r, Lu_plus, Lu_minus, U0, V0, S0, alpha, gamma);
    V0 = updateU(R', UVT',w, r, Lv_plus, Lv_minus, V0, U0, S0', alpha, lambda);
    S0 = updateS(S0, R, UVT, w, r, U0, V0, alpha);
    ite = ite + 1;
end

U = U0;
V = V0;
S = S0;
end

function [UVT] = get_UVT(R, U, S, V)

[m, n] = size(R);

[I, J] = find(R);
iSize = size(I, 1);
K = ones(iSize,1);
count = 0;
US=U*S;
for j = 1:n
    id = find(R(:,j));
    len = length(id);
    K(count+1:count+len) = US(id, :) * V(j,:)';
    count = count + len;
end

UVT = sparse(I, J, K, m, n);

end

function [U1] = updateU(R, UVT, w, p, Lu_plus, Lu_minus, U0, V, S0, lambda, gamma)

[m,n]=size(R);
U1 = U0 .* sqrt( ((1-w*p)*R*V*S0' + ones(m,1)*p*((w*ones(1,n))*V*S0') + gamma .* Lu_minus * U0) ./ ((1-w)*UVT*V*S0' + w * (U0*S0*(V'*V)*S0')  + gamma.* Lu_plus * U0 + lambda * U0) );

U1(isnan(U1)) = 0;

end
function [S1] = updateS(S0, R, UVT, w, p, U, V,lambda)
[m,n]=size(R);
S1 = S0.*sqrt( ((1-w*p)*(U')*R*V + (U')*(w*p*ones(m,n))*V )./( (1-w)*(U')*UVT*V + w*(U')*U*S0*(V')*V+lambda*S0) );
S1(isnan(S1)) = 0;
end
