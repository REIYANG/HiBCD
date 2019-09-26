function [S] = S_init_FPG(S,A,F,G,Y_H,Y_M)
[M,N] = size(A);
uk = 1;
S_prev = S;
AtA = A'*A;
AtY_H = A'*Y_H;
FA = F*A;
[U,~,~] = svd(eye(N)-1/N*ones(N,1)*ones(1,N));
Phi = U(:,1:N-1);
APhi = A*Phi;
thetaG_FFt = PowerMethod(G'*G,100)*eye(M)+F'*F;
delta = eps;
beta_S = max(norm(APhi'*thetaG_FFt*APhi,2),delta);
for i = 1:100
    uk_next = (1+sqrt(1+4*uk^2))/2;
    iota_k = (uk-1)/uk_next;
    uk = uk_next;
    % -------- S-subproblem
    S_ex = S+iota_k*(S-S_prev);
    S_prev = S;
    nabla_S = (AtA*(S_ex*G)-AtY_H)*G'+FA'*(FA*S_ex-Y_M);
    
    
    S = Simplex_Proj(S_ex-1/beta_S*nabla_S);
end
end

function maxEig = PowerMethod(M,iteraNum)
n = size(M,1); x = rand(n,1)-0.5; x = x./abs(x);
for k = 1 : iteraNum ;
    x = M*x;
end
maxEig = (x'*M*x)/(x'*x);
end