function TriFacREMAP_timecomplexity()
%[1] Hansaim Lim, Aleksandar Poleksic, Hanghang Tong, Yuan Yao, Di He, Luke Zhuang, Patrick Meng, and Lei Xie, "Large-scale Off-target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing" , under review
%chem_chem_mat and prot_prot_mat must be square matrices
%number of rows and number of columns of R must match the number of chemicals and proteins, respectively
%prediction matrix Y does NOT show scores for known associations
%in other words, scores for known association pairs are set to 0
%Please refer to the paper for detail
maxNumCompThreads(2); %determine the maximum number of cores to use
%load('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\thesis work\TriFacREMAP\script\testpack\datamat\ChEMBLCYP450\chem_chem.mat')
%load('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\thesis work\TriFacREMAP\script\testpack\datamat\ChEMBLCYP450\chem_prot.mat')
%load('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\thesis work\TriFacREMAP\script\testpack\datamat\ChEMBLCYP450\prot_prot.mat')
% load('.\testpack\datamat\kinase\chem_prot.mat')
% load('.\testpack\datamat\kinase\chem_chem.mat')
% load('.\testpack\datamat\kinase\prot_prot.mat')


%for ZINC dataset
load('..\..\external data\ZINC\chem_prot_zinc.mat')
load('..\..\external data\ZINC\chem_chem_zinc.mat')
load('..\..\external data\ZINC\prot_prot_zinc.mat')
chem_prot=sparse(chem_prot_zinc);
chem_chem=sparse(chem_chem_zinc);
prot_prot=sparse(prot_prot_zinc);
%for ZINC dataset


chem_prot=chem_prot(sum(chem_prot,2)>0,:);
chem_chem=chem_chem(sum(chem_prot,2)>0,sum(chem_prot,2)>0);
chem_prot=chem_prot(:,sum(chem_prot,1)>0); %size=20487,27
prot_prot=prot_prot(sum(chem_prot,1)>0,sum(chem_prot,1)>0);
nchem=size(chem_prot,1);
nprot=size(chem_prot,2);

K=4; %repeat each condition n times to get stdev
outfile='./TREMAP_Timecomplexity_stdev_1090I_large_new.txt';
fid=fopen(outfile,'at+');
fprintf(fid,'%s\t%3d\n','N=',K);
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
    'numDrug','numTarget','lowDrug','lowTarget','a*b*c*d','Avg. time(s)','Std. time(s)');

	% para: p_reg, squared p_weight, p_imp, rank_u, rank_v, p_iter, p_chem, p_prot

% % range1=2:1:5; %column range
% % range2=2:1:5; %row range
range2=200000:200000:1000000;
range1=100000:20000:200000;
rank_u=400;
rank_v=200;
for i1=1:numel(range1)
    r1=range1(i1);
    colsample=datasample(1:nprot,r1);
    R1=chem_prot(:,colsample);
    N=prot_prot(colsample,colsample);
    for i2=1:numel(range2)
        r2=range2(i2);
        rowsample=datasample(1:nchem,r2);
        R=R1(rowsample,:);
        M=chem_chem(rowsample,rowsample);
        %get number of chemical and protein
        Times=zeros(1,K);
        for k=1:K
            tic;
            m=size(M, 1);
            n=size(N, 1);

            summ = sum(M,2); %sum by rows
            Dm = spdiags(summ,0,m,m);
            Lu = Dm - M;

            sumn = sum(N,2); %sum by rows
            Dn = spdiags(sumn,0,n,n);
            Lv = Dn - N;
            para = [0.1, 0.1, 0.01, rank_u, rank_v, 100, 0.5, 0.5];
            [U, S, V] = updateUSV(R, Lu, Lv, para);
            t=toc;
            Times(1,k)=t;
            clear U S V Dm Dn Lu Lv t
        end
        Avg_t=mean(Times);
        Std_t=std(Times);
        fprintf(fid,'%5d\t%5d\t%5d\t%5d\t%15d\t%.8g\t%.8g\n',...
            r2,r1,rank_u,rank_v, r2*r1*rank_u*rank_v ,Avg_t,Std_t);
    end
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
