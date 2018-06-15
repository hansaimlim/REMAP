function Figure3_example()
%[1] Lim, H., Poleksic, A., Yao, Y., Tong, H., He, D., Zhuang, L., Meng, P. and Xie, L., 2016. Large-Scale Off-Target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing. PLOS Comput Biol, 12(10), p.e1005135.
maxNumCompThreads(2) %determine the maximum number of cores to use
para = [0.1, 0.1, 0.01, 200, 300, 0.75, 0.1];	% para: p_reg, squared p_weight, p_imp, rank, p_iter, p_chem, p_prot
kFold=10;
output_file='./REMAP_fig3_example.txt';
fileid=fopen(output_file,'a+t');
%chem_chem_zinc and protein_protein_zinc_blast matrices from chem-chem and prot-prot files
load ../benchmark/sim/chem_chem_zinc;
load ../benchmark/sim/prot_prot_zinc;

%cp3=csvread('../benchmark/NTNL/train_N2L16to20.csv');
%cp4=csvread('../benchmark/NTNL/test_N2L16to20.csv');
%get number of chemical and protein
m=size(chem_chem_zinc, 1);
n=size(prot_prot_zinc, 1);
%convert csv to matrix


summ = sum(chem_chem_zinc,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_zinc;

sumn = sum(prot_prot_zinc,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_zinc;

trainfiles=["../benchmark/NTNL/train_N2L1to5.csv","../benchmark/NTNL/train_N2L6to10.csv",...
    "../benchmark/NTNL/train_N2L11to15.csv","../benchmark/NTNL/train_N2L16to20.csv",...
    "../benchmark/NTNL/train_N2L21more.csv"];
testfiles=["../benchmark/NTNL/test_N2L1to5.csv","../benchmark/NTNL/test_N2L6to10.csv",...
    "../benchmark/NTNL/test_N2L11to15.csv","../benchmark/NTNL/test_N2L16to20.csv",...
    "../benchmark/NTNL/test_N2L21more.csv"];

for fidx=1:length(trainfiles)
trainfile=trainfiles(fidx);
testfile=testfiles(fidx);
cptrain=csvread(trainfile);
cptest=csvread(testfile);
train = sparse(cptrain(:,1),cptrain(:,2),1,m,n);%12384 chemicals and 3500 proteins in ZINC
cv=cvpartition(length(cptest),'Kfold',kFold);

TPR=zeros(kFold,1);
    for k=1:kFold
      tic;
      %add kfold test pairs into train
      trainpairs=cptest(training(cv,k),:);
      trainmat=sparse(trainpairs(:,1),trainpairs(:,2),1,m,n);
      testpairs=cptest(test(cv,k),:);
      train_k=train+trainmat;
      test_k=sparse(testpairs(:,1),testpairs(:,2),1,m,n);
      %mask kfold test pairs in test
      %evaluate on test matrix
      [U, V] = updateUV(train_k, Lu, Lv, para);
      P=U*V';P=P(sum(test_k,2)>0,:);
      testcompressed=test_k(sum(test_k,2)>0,:);
      tprbyrow = TPRbyRowRank(FindTrues(P, testcompressed), 40);   %max cutoff rank 100
      tpr=tprbyrow(35,2);
      TPR(k,1)=tpr;
      disp(['Rank=' num2str(para(4)) ' Iter=' num2str(para(5)) ' TPR at top 35th=' num2str(tpr) ])
      toc
      clear P U V testcompressed
    end
text=['file=' testfile ': Top35 avgTPR=' num2str(mean(TPR)) ' stdTPR=' num2str(std(TPR)) '\n'];
fwrite(fileid, strjoin(text));
end

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
