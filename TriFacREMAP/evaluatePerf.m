function [MAP,HLU,MPR]=evaluatePerf(test,predictmat)

m = size(predictmat, 1);
n = size(predictmat, 2);

MAP = zeros(m,1);
HLU_1 = zeros(m,1);
HLU_2 = zeros(m,1);
PR_1 = zeros(m,1);
PR_2 = zeros(m,1);

for i = 1:m
    real = test(i,:)';
    predict = predictmat(i,:)';
    
    MAP(i,1) = computeAP(real, predict);
    [HLU_1(i,1), HLU_2(i,1)] = computeHLU(real, predict);
    [PR_1(i,1), PR_2(i,1)] = computePR(real, predict);
end

[I] = find(MAP>0);
MAP = MAP(I);
MAP = mean(MAP);
HLU = 100 * sum(HLU_1) / sum(HLU_2);
MPR = sum(PR_1) / sum(PR_2);
end