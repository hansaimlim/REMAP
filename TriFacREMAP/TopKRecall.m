function [topKrec] = TopKRecall(PredMat, TrueMat, maxK)
%get average TPR by cutoff Row rank (top K)
%output contains vector of dimension (maxK, 3)

if size(PredMat,1)~=size(TrueMat,1)||size(PredMat,2)~=size(TrueMat,2)
    msg='The PredMat and TrueMat must be in same size.';
    error(msg)
end
if maxK > size(PredMat,2)
   maxK=size(PredMat,2); %number of predictions smaller than K 
end
n=size(PredMat,1);
TopKRecRows=zeros(n,maxK); %top K recall for each row
topKrec=zeros(3,maxK); %[k, avgRec, stdRec]
%TopKPrecRows=zeros(n,maxK);
for i=1:n
    [~,order]=sort(PredMat(i,:),'descend');    
    TrueIdx=find(TrueMat(i,:));
    cp=numel(TrueIdx); %condition positive
    for k=1:maxK
         if ismember(order(k),TrueIdx)
             %add true count for the k
             if k==1
                 TopKRecRows(i,k)=1; %the top rank
             else
                 TopKRecRows(i,k)=TopKRecRows(i,k-1)+1;
             end
         else
             if k==1
                 TopKRecRows(i,k)=0;
             else
                 TopKRecRows(i,k)=TopKRecRows(i,k-1);
             end
         end
    end
    TopKRecRows(i,:)=TopKRecRows(i,:)/cp;
end
topKrec(1,:)=1:maxK;
topKrec(2,:)=mean(TopKRecRows);
topKrec(3,:)=std(TopKRecRows);
end
