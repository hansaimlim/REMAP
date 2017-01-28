function P=COSINE_Optimization(chem_prot,chem_chem,prot_prot)
% Lim, H., Gray, P., Xie, L., & Poleksic, A. (2016). 
% Improved genome-scale multi-target virtual screening via a novel collaborative filtering approach to cold-start problem. 
% Scientific Reports, 6, 38860. http://doi.org/10.1038/srep38860
maxNumCompThreads(3); %determine the maximum number of cores to use

R = chem_prot;
M = chem_chem;
N = prot_prot;

kFold=10;
fid=fopen('../COSINE_optimization_ZINC.txt','at+');
ColdStartRowIdx=find(sum(R,2)==1);
%    ColdStartColIdx=find(sum(R,1)==1);

ranks=100:100:500;
lRs=0:0.33:1.0;
lMs=0:0.33:1.0;
lNs=0:0.33:1.0;
total_grid=numel(lRs)*numel(lMs)*numel(lNs)+numel(ranks);
grid_searched=0;

J=3;
[DM, nM]= GetDiag(M,J);
[DN, nN]= GetDiag(N,J);

DMM = DM-nM;
DNN = DN-nN;  

cv=cvpartition(length(ColdStartRowIdx),'Kfold',kFold);
best_lambdas=zeros(1,3); %[lR, lM, lN]
best_auc=0;
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n','rank','lR','lM','lN','avgAUC','stdAUC');    
for ilr=1:numel(lRs)
    for ilm=1:numel(lMs)
        for iln=1:numel(lNs)
            rnk=200;wp=0.1;lR=lRs(ilr);lM=lMs(ilm);lN=lNs(iln);
            paras=[rnk,wp,lR,lM,lN];
            tic;
            [avgauc,stdauc]=COSINE_10CV(ColdStartRowIdx,cv,kFold,R,paras,DMM,DNN,nM,nN);
            fprintf(fid,'%5d\t%.3g\t%.3g\t%.3g\t%.5g\t%.5g\n',rnk,lR,lM,lN,avgauc,stdauc);
            if avgauc>best_auc
                best_auc=avgauc;
                best_lambdas=[lR,lM,lN];        
            end
            grid_searched=grid_searched+1;
            disp(['COSINE_Optimization Grid Search - ' num2str(grid_searched) ...
                ' out of ' num2str(total_grid) ' grids searched.']);
            toc
        end
    end
end

best_paras=[200,best_lambdas(1),best_lambdas(2),best_lambdas(3)];
for ir=1:numel(ranks)
        rnk=ranks(ir);lR=best_lambdas(1);lM=best_lambdas(2);lN=best_lambdas(3);
        paras=[rnk,lR,lM,lN];
        tic;
        [avgauc,stdauc]=COSINE_10CV(ColdStartRowIdx,cv,kfold,R,paras,DMM,DNN,nM,nN);
        fprintf(fid,'%5d\t%.3g\t%.3g\t%.3g\t%.5g\t%.5g\n',rnk,lR,lM,lN,avgauc,stdauc);
        if avgauc>best_auc
            best_auc=avgauc;
            best_paras(1)=rnk;      
        end
        grid_searched=grid_searched+1;
        disp(['COSINE_Optimization Grid Search - ' num2str(grid_searched) ...
            ' out of ' num2str(total_grid) ' grids searched.']);
        toc
end
fclose(fid);
disp(['Best avg. AUC=' num2str(best_auc) ' at rank=' num2str(best_paras(1)) ...
    ', lR=' num2str(best_paras(2)) ', lM=' num2str(best_paras(3)) ...
    ', lN=' num2str(best_paras(4))]);

%get prediction matrix P using the best parameters
W = max(1, 6 * Train);
Q = zeros(m,n);  
P=COSINE(R,DMM,DNN,nM,nN,W,Q,best_paras(2),best_paras(3),best_paras(4),best_paras(1));
    
end

function [avgauc,stdauc]=COSINE_10CV(ColdStartRowIdx,cv,kFold,R,paras,DMM,DNN,nM,nN)
    AUC=zeros(1,kFold);
    rnk=paras(1);lR=paras(2);lM=paras(3);lN=paras(4);
    [m,n] = size(R);
    for k=1:kFold
        %Cross validation, K-fold
        testidx=ColdStartRowIdx(test(cv,k)); %k-th test set
        Train=R;
        Train(testidx,:)=0; %mask test pairs to 0
        Test=R-Train;

        W = max(1, 6 * Train);
        Q = zeros(m,n);            

        EXC=COSINE(Train,DMM,DNN,nM,nN,W,Q,lR,lM,lN,rnk);
        [~,~,~,auc]=perfcurve(Test(:),EXC(:),1);
        AUC(1,k)=auc;
    end
    avgauc=mean(AUC);
    stdauc=std(AUC);
end

function EXC=COSINE(Train,DMM,DNN,nM,nN,W,Q,lR,lM,lN,rnk)
iter=200;
TANIM = 0.7; FGSCORE = 0.95; WGHT = 7; WP=1.0;
J = 3;
QUICK = 1;
M_cut = 0.3;
N_cut = 0.5;
[ExcludedRows, unimportant] = find(sum(Train,2)==0); 
[unimportant, ExcludedColumns] = find(sum(Train,1)==0); 
[F, G] = WeightImputeLogFactorization(Train,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);
[nF, nG, HI_IND] = WeightedProfile(F, G, nM, nN, ExcludedRows, ExcludedColumns, J + 2, TANIM, WP, M_cut, N_cut);

    EXC = GetP(nF*nG');
    if QUICK == 0
        LOW = EXC < FGSCORE;
        HIGH = EXC >= FGSCORE;
        EXC(HIGH) = 1;
        EXC(LOW) = 0;
        W(HI_IND > 0, :) = WGHT;
        EXC = max(EXC,Train);

        [F, G] =  WeightImputeLogFactorization(EXC,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);  

        [nF, nG, HI_IND] = WeightedProfile(F, G, nM, nN, ExcludedRows, ExcludedColumns, J + 2, TANIM, WP, M_cut, N_cut);
        EXC = GetP(nF*nG'); 
    end
end
