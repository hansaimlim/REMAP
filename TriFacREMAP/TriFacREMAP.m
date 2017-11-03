function [U,S,V]=TriFacREMAP(R, chem_chem_mat, prot_prot_mat,param)
%obtain low-rank matrices using optimized parameters
%maxNumCompThreads(2); %determine the maximum number of cores to use
default_reg=0.1;
default_weight=0.1;
default_impute=0.1;
default_iter=100;
if nargin < 4
    para = [default_reg, default_weight, default_impute, 200, 15, ...
        default_iter, 0.2, 0.5];	% para: p_reg, squared p_weight, p_imp, rank_u, rank_v, p_iter, p_chem, p_prot
end
if length(param) == 4
    para(4)=param(1);para(5)=param(2);para(7)=param(3);para(8)=param(4);
elseif length(param) == 8
    para=param;
else
    msg=['The parameter must contain 4 or 8 values.\n' ...
        'e.g.)para=[0.1, 0.1, 0.01, 200, 50, 300, 0.75, 0.1]\n'...
        'or para=[200, 50, 0.75, 0.25] (rank_u, rank_v, p_chem, p_prot)\n'];
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
[U, S, V] = updateUSV(R, Lu, Lv, para);

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
