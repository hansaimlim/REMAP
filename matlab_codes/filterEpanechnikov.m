function Mat=filterEpanechnikov(P,R,probcutoff)
%Find predicted associations from P with high score
%create Epanechnikov kernel for each row
%rows are independent

[m,n]=size(P);
Mat=zeros(floor(m*n/5),3); %[rowIdx,colIdx,cdf]
count=1;
Epa=zeros(m,n);
for i=1:m
   row=P(i,:)';
   known=R(i,:)';
   rownonzero=row(~known);
   pd=fitdist(rownonzero,'Kernel','Kernel','Epanechnikov');
   cd=cdf(pd,row);
   Epa(i,:)=cd;
   cols=find(cd>=probcutoff);
   for j=1:numel(cols)
       colIdx=cols(j);
       if known(colIdx)==1
           continue
       end
       prob=cd(colIdx);
       Mat(count,:)=[i,colIdx,prob];
       count=count+1;
   end
end
Mat=Mat(1:count-1,:);
end