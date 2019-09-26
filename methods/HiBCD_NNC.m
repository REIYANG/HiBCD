function [Y_HiBCD,f_record,time_record,count] = HiBCD_nuclear(HSI,MSI,F,G,N,varargin)
S_init = [];
A_init = [];
stopThreshd = 1e-5;
maxIter = 500;
gamma = 20;
% --------read the optional parameters
if (nargin-length(varargin))~=5
    error('Wrong number of required parameters');
end
if rem(length(varargin),2)==1
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'MAXITER'
                maxIter = varargin{i+1};
            case 'STOPPING'
                stopThreshd = varargin{i+1};
            case 'GAMMA'
                gamma = varargin{i+1};
            case 'S_INIT'
                S_init = varargin{i+1};
            case 'A_INIT'
                A_init = varargin{i+1};
            case 'A_UPDATE'
                A_flag = varargin{i+1};
            case 'S_UPDATE'
                S_flag = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end
% -------- load parameters
M_M = size(MSI,3);
M = size(HSI,3);
Y_M = reshape(MSI,[],M_M)';
Y_H = reshape(HSI,[],M)';
[L1,L2,~] = size(MSI);
delta = eps;
% -------- frequently used var
theta_F = norm(F*F',2);
if S_flag ~= 2
    Phi = [eye(N-1);-ones(1,N-1)];
    thetaG_FFt = PowerMethod(G'*G,100)*eye(M)+F'*F;
end
% -------- algorithm initialization
if isempty(A_init)||isempty(S_init)
    A_init = max(min(SPA(Y_H,N),1),0);
    S_init = Simplex_Proj(pinv(F*A_init)*Y_M);
end
for n = 1:N
    [U,Sig,V] = svds(reshape(S_init(n,:),L1,L2),L1);
    diag_Sig = diag(Sig);
    if sum(diag_Sig) <= gamma
    else
        diag_Sig = gamma*Simplex_Proj(diag_Sig/gamma);
    end
    S_n = U(:,1:length(diag_Sig))*diag(diag_Sig)*V(:,1:length(diag_Sig))';
    S_init(n,:) = S_n(:)';
end
S = S_init;
S_prev = S;
A = max(min(A_init,1),0); A_prev = A;
uk = 0;
time_record = zeros(maxIter+1,1);
f_record = zeros(maxIter+1,1);
tic;
FA = F*A;
YM_FAS = Y_M-FA*S;
SG = S*G;
YH_ASG = Y_H-A*SG;
f_curr = 1/2*(YH_ASG(:)'*YH_ASG(:)+YM_FAS(:)'*YM_FAS(:));
f_record(1) = f_curr;
U = rand(L1,N);
V = rand(L2,N);
% -------- HiBCD
for count = 1:maxIter
    if S_flag==3||A_flag==3
        uk_next = (1+sqrt(1+4*uk^2))/2;
        iota_k = (uk-1)/uk_next;
        uk = uk_next;
    end
    % -------- S-subproblem
    AtA = A'*A;
    AtY_H = A'*Y_H;
    if S_flag == 3 % with extrapolation
        S_ex = S+iota_k*(S-S_prev);
        S_prev = S;
    else % without extrapolation
        S_ex = S;
    end
    nabla_S = (AtA*(S_ex*G)-AtY_H)*G'+FA'*(FA*S_ex-Y_M);
    if S_flag == 2
        P = zeros(size(S_ex));
        for n = 1:N
            nabla_Sn = reshape(nabla_S(n,:),L1,L2);
            [U(:,n),V(:,n)] = PowerSVD(nabla_Sn,5,V(:,n));
            vec_tmp = gamma*U(:,n)*V(:,n)';
            P(n,:) = -vec_tmp(:)';
        end
        D_S = P-S_ex;
        ADG = A*(D_S*G); FAD = FA*D_S;
        alpha_k = min(1,-(nabla_S(:)'*D_S(:))/...
            (ADG(:)'*ADG(:)+FAD(:)'*FAD(:)+delta*D_S(:)'*D_S(:)));
        S = S_ex+alpha_k*D_S;
    else
        beta_S = max(norm(A'*thetaG_FFt*A,2),delta);
        S = NuclearProj(S_ex-1/beta_S*nabla_S,L1,L2,gamma);
    end
    % -------- A-subproblem
    SG = S*G ;
    if A_flag == 3 % with extrapolation
        A_ex = A+iota_k*(A-A_prev);
        A_prev = A;
    else % without extrapolation
        A_ex = A;
    end
    nabla_A = (A_ex*SG-Y_H)*SG'+F'*(((F*A_ex)*S-Y_M)*S');
    if A_flag == 2 % FW update
        D_A = (nabla_A<0)-A_ex;
        DSG = D_A*SG; FDS = (F*D_A)*S;
        alpha_k = min(1,-(nabla_A(:)'*D_A(:))/...
            (DSG(:)'*DSG(:)+FDS(:)'*FDS(:)+delta*D_A(:)'*D_A(:)));
        A = A_ex+alpha_k*D_A;
    else % (F)PG update
        beta_A = max(norm(theta_F*(S*S')+SG*SG',2),delta);
        A = min(1,max(0,A_ex-1/beta_A*nabla_A));
%         if sum(sum(nabla_A.*(A-A_prev)))>0
%             uk = 1;
%         end
    end
    
    % -------- check convergence
    f_prev = f_curr;
    FA = F*A;
    SG = S*G;
    YM_FAS = Y_M-FA*S;
    YH_ASG = Y_H-A*SG;
    f_curr = 1/2*(YH_ASG(:)'*YH_ASG(:)+YM_FAS(:)'*YM_FAS(:));
    time_record(count+1) = toc+time_record(count);
    tic;
    f_record(count+1) = f_curr;
    obj_chng = abs(f_prev-f_curr)/f_prev;
    if count>=maxIter||obj_chng<stopThreshd
        break ;
    end
end
f_record(count+2:end) = [];
time_record(count+2:end) = [];
f_record(end);
for n = 1:N
    [U,Sig,V] = svds(reshape(S(n,:),L1,L2),L1);
    diag_Sig = diag(Sig);
    if sum(diag_Sig) <= gamma+1e-10
    else
        n
    end
end
Y_HiBCD = A*S;
end

function [u,v] = PowerSVD(M,iteraNum,x)
MM = M'*M;
for k = 1:iteraNum
    x = MM*x;
end
v = x/norm(x);
Mv = M*v;
u = Mv/norm(Mv);
end

function S_proj = NuclearProj(S,L1,L2,gamma)
[N,L] = size(S);
S_proj = zeros(N,L);
for n = 1:N
    Sn = reshape(S(n,:),L1,L2);
    [U,Sig,V] = svds(Sn,N);
    diag_Sig = diag(Sig);
    if sum(diag_Sig)<=gamma
        Sn_proj = Sn;
    else
        Sn_proj = U(:,1:length(diag_Sig))*(diag(Simplex_Proj(diag_Sig/gamma)*gamma))*V(:,1:length(diag_Sig))';
    end
    S_proj(n,:) = Sn_proj(:)';
end
end

function maxEig = PowerMethod(M,iteraNum)
n = size(M,1); x = rand(n,1)-0.5; x = x./abs(x);
for k = 1 : iteraNum 
    x = M*x;
    tmp = x'*x;
    if tmp>=1e200||tmp<=1e-200
        x = x/norm(x);
    end
end
maxEig = (x'*M*x)/(x'*x);
end

function S_Proj = Simplex_Proj(S)
S_Proj = max(bsxfun(@minus,S,max(bsxfun(@rdivide,cumsum(sort(S,1,'descend'),1)-1,(1:size(S,1))'),[],1)),0);
end