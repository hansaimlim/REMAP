function optParams=TriFacREMAP_optimization(R, chem_chem_mat, prot_prot_mat,outfile)
%Optimize parameters by cross-validation
%maxNumCompThreads(2); %determine the maximum number of cores to use
kFold = 10;
fid=fopen(outfile,'at+');
optParams=[0 0 0 0]; %rank_u, rank_v, p6, p7
optTopKRecall=0;
%get number of chemical and protein
m=size(chem_chem_mat, 1);
n=size(prot_prot_mat, 1);

summ = sum(chem_chem_mat,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_mat;

sumn = sum(prot_prot_mat,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_mat;
maxrank=max(m,n);
minrank=min(m,n);

if maxrank>1500
    ranks_u=floor(maxrank/100):floor(maxrank/70):floor(maxrank/10);
else
    ranks_u=floor(maxrank/10):floor(maxrank/20):floor(maxrank/3);
end

ranks_v=floor(minrank/10):floor(minrank/10):floor(minrank/3);
% iters=100:100:200;
p7s=0:0.33:1.0;
p8s=0:0.33:1.0;


%prepare CV object
positive_idx=find(R>0);
if n<=5
    maxK=5;
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'LowRank_U','LowRank_V','p7','p8','Top1avgRec','Top1stdRec','Top2avgRec',...
        'Top2stdRec','Top3avgRec','Top3stdRec','Top4avgRec','Top4stdRec');    
elseif n<=30
    maxK=20;
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'LowRank_U','LowRank_V','p7','p8','Top5avgRec','Top5stdRec','Top10avgRec',...
        'Top10stdRec','Top15avgRec','Top15stdRec','Top20avgRec','Top20stdRec');
elseif n<=50
    maxK=30;
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'LowRank_U','LowRank_V','p7','p8','Top5avgRec','Top5stdRec',...
        'Top10avgRec','Top10stdRec','Top15avgRec','Top15stdRec',...
        'Top25avgRec','Top25stdRec');
elseif n<=100
    maxK=50;
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'LowRank_U','LowRank_V','p7','p8','Top5avgRec','Top5stdRec',...
        'Top10avgRec','Top10stdRec','Top30avgRec','Top30stdRec',...
        'Top40avgRec','Top40stdRec');    
elseif n<=150
    maxK=100;
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'LowRank_U','LowRank_V','p7','p8','Top5avgRec','Top5stdRec',...
        'Top10avgRec','Top10stdRec','Top50avgRec','Top50stdRec',...
        'Top70avgRec','Top70stdRec'); 
elseif n<=200
    maxK=150;
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'LowRank_U','LowRank_V','p7','p8','Top10avgRec','Top10stdRec',...
        'Top30avgRec','Top30stdRec','Top70avgRec','Top70stdRec',...
        'Top120avgRec','Top120stdRec'); 
else
    maxK=200;
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'LowRank_U','LowRank_V','p7','p8','Top10avgRec','Top10stdRec',...
        'Top40avgRec','Top40stdRec','Top100avgRec','Top100stdRec',...
        'Top150avgRec','Top150stdRec'); 
end
cv=cvpartition(length(positive_idx),'Kfold',kFold);

for i1=1:numel(ranks_u)
    for i2=1:numel(ranks_v)
        for i3=1:numel(p7s)
            for i4=1:numel(p8s)
                rank_u=ranks_u(i1);
                rank_v=ranks_v(i2);
                p6=p7s(i3); %p_chem; written as p6 for consistency with REMAP
                p7=p8s(i4); %p_prot; written as p7 for consistency with REMAP
                para = [0.1, 0.1, 0.01, rank_u, rank_v, 100, p6, p7];
                topKrec=zeros(kFold,maxK);
                for k=1:kFold
                    %Cross validation, K-fold
                    testidx=positive_idx(test(cv,k)); %k-th test set
                    Train=R;
                    Train(testidx)=0; %mask test pairs to 0
                    Test=R-Train;
                    TestCompressed=Test(sum(Test,2)>0,:);%only rows with nonzero element
                    [U, S, V] = updateUSV(Train, Lu, Lv, para);
                    Pred=U*S*V';
                    Pred=Pred(sum(Test,2)>0,:); %same rows as test compressed matrix
                    topkrec=TopKRecall(Pred,TestCompressed,350);
                    topKrec(k,1:size(topkrec,2))=topkrec(2,:);
                    clear U S V Pred
                end
                avgTopKRec=mean(topKrec);
                stdTopKRec=std(topKrec);
                if n<=5
                    fprintf(fid,'%5d\t%5d\t%.3g\t%.3g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n',...
                        ranks_u(i1),ranks_v(i2),p7s(i3),p8s(i4),avgTopKRec(1),stdTopKRec(1),avgTopKRec(2),stdTopKRec(2),...
                        avgTopKRec(3),stdTopKRec(3),avgTopKRec(4),stdTopKRec(4));
                    selectedTopKRec=avgTopKRec(2);
                elseif n<=30
                    fprintf(fid,'%5d\t%5d\t%.3g\t%.3g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n',...
                        ranks_u(i1),ranks_v(i2),p7s(i3),p8s(i4),avgTopKRec(5),stdTopKRec(5),avgTopKRec(10),stdTopKRec(10),...
                        avgTopKRec(15),stdTopKRec(15),avgTopKRec(20),stdTopKRec(20));
                    selectedTopKRec=avgTopKRec(10);
                elseif n<=50
                    fprintf(fid,'%5d\t%5d\t%.3g\t%.3g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n',...
                        ranks_u(i1),ranks_v(i2),p7s(i3),p8s(i4),avgTopKRec(5),stdTopKRec(5),avgTopKRec(10),stdTopKRec(10),...
                        avgTopKRec(15),stdTopKRec(15),avgTopKRec(25),stdTopKRec(25));
                    selectedTopKRec=avgTopKRec(15);
                elseif n<=100
                    fprintf(fid,'%5d\t%5d\t%.3g\t%.3g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n',...
                        ranks_u(i1),ranks_v(i2),p7s(i3),p8s(i4),avgTopKRec(5),stdTopKRec(5),avgTopKRec(10),stdTopKRec(10),...
                        avgTopKRec(30),stdTopKRec(30),avgTopKRec(40),stdTopKRec(40));
                    selectedTopKRec=avgTopKRec(30);
                elseif n<=150
                    fprintf(fid,'%5d\t%5d\t%.3g\t%.3g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n',...
                        ranks_u(i1),ranks_v(i2),p7s(i3),p8s(i4),avgTopKRec(5),stdTopKRec(5),avgTopKRec(10),stdTopKRec(10),...
                        avgTopKRec(50),stdTopKRec(50),avgTopKRec(70),stdTopKRec(70));
                    selectedTopKRec=avgTopKRec(50);
                elseif n<=200
                    fprintf(fid,'%5d\t%5d\t%.3g\t%.3g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n',...
                        ranks_u(i1),ranks_v(i2),p7s(i3),p8s(i4),avgTopKRec(10),stdTopKRec(10),avgTopKRec(30),stdTopKRec(30),...
                        avgTopKRec(70),stdTopKRec(70),avgTopKRec(120),stdTopKRec(120));
                    selectedTopKRec=avgTopKRec(70);
                else
                    fprintf(fid,'%5d\t%5d\t%.3g\t%.3g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n',...
                        ranks_u(i1),ranks_v(i2),p7s(i3),p8s(i4),avgTopKRec(10),stdTopKRec(10),avgTopKRec(40),stdTopKRec(40),...
                        avgTopKRec(100),stdTopKRec(100),avgTopKRec(150),stdTopKRec(150));
                    selectedTopKRec=avgTopKRec(100);
                end
                if selectedTopKRec>optTopKRecall
                    optParams=[rank_u, rank_v, p6, p7]; %current best rank, p6, p7
                end
            end
        end
    end
end
             
fprintf(fid,'%s\t%s\t%s\t%s\n','OptRank_u','OptRank_v','Optp6','Optp7');
fprintf(fid,'%5d\t%5d\t%.3g\t%.3g\n',optParams(1),optParams(2),optParams(3),optParams(4));
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
