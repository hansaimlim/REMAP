function get_train_mat()
%mask test and validation dataset
cd('D:\OneDrive - Cuny GradCenter\thesis work\TriFacREMAP\script\testpack\datamat\PubChemCYP450\')
load('chem_prot.mat')
testPyidx=csvread('chem_prot_test_pyIndex.csv');
validPyidx=csvread('chem_prot_valid_pyIndex.csv');
chem_prot_train=chem_prot;
for i=1:size(testPyidx,1)
   chemidx=testPyidx(i,1)+1; %+1 for python indices
   protidx=testPyidx(i,2)+1;
   chem_prot_train(chemidx,protidx)=0; %mask
end

for i=1:size(validPyidx,1)
   chemidx=validPyidx(i,1)+1; %+1 for python indices
   protidx=validPyidx(i,2)+1;
   chem_prot_train(chemidx,protidx)=0; %mask
end
save('chem_prot_train.mat','chem_prot_train');

end