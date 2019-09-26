function [Y_HiBCD,f_record,time_record,count] = HiBCD(HSI,MSI,F,G,N,varargin)
S_init = [];
A_init = [];
stopThreshd = 1e-5;
maxIter = 500;
tight_flag = 1;
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
            case 'S_INIT'
                S_init = varargin{i+1};
            case 'A_INIT'
                A_init = varargin{i+1};
            case 'A_UPDATE'
                A_flag = varargin{i+1};
            case 'S_UPDATE'
                S_flag = varargin{i+1};
            case 'TIGHT'
                tight_flag = varargin{i+1};
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
[L_1,L_2,~] = size(MSI);
L = L_1*L_2;
delta = eps;
% -------- frequently used var 
thetaG_FFt = PowerMethod(G'*G,100)*eye(M)+F'*F;
theta_F = norm(F*F',2); 
% Phi = [eye(N-1);-ones(1,N-1)];
[U,~,~] = svds(eye(N)-1/N*ones(N,1)*ones(1,N),N);
Phi = U(:,1:N-1);
% -------- algorithm initialization 
if isempty(A_init)
    A_init = init_AS(Y_H,N);
end
if isempty(S_init)
    S_init = rand(N,L);
    S_init = bsxfun(@rdivide,S_init,sum(S_init));
end
S = S_init; S_prev = S;
A = A_init; A_prev = A;
uk = 1; 
time_record = zeros(maxIter+1,1);
f_record = zeros(maxIter+1,1);
tic;
FA = F*A;
YM_FAS = Y_M-FA*S;
SG = S*G;
YH_ASG = Y_H-A*SG;
f_curr = 1/2*(YH_ASG(:)'*YH_ASG(:)+YM_FAS(:)'*YM_FAS(:));
f_record(1) = f_curr;
% -------- HiBCD
for count = 1:maxIter
    if S_flag==3||A_flag==3
        uk_next = (1+sqrt(1+4*uk^2))/2;
        iota_k = (uk-1)/uk_next;
        uk = uk_next;
    end
    % -------- S-subproblem
    if S_flag == 3 % with extrapolation
        S_ex = S+iota_k*(S-S_prev);
        S_prev = S;
    else % without extrapolation
        S_ex = S;
    end
    AtA = A'*A;
    AtY_H = A'*Y_H;
    if S_flag == 2 % FW update
        nabla_S = (AtA*(S_ex*G)-AtY_H)*G'+FA'*(FA*S_ex-Y_M);
        [~,colIndex] = min(nabla_S);
%         D_S = ind2vec(colIndex,N)-S_ex;
        tmp = sparse(1:numel(colIndex), colIndex,1)';
        D_S = [tmp; sparse(N-size(tmp,1),L)]-S_ex;
        ADG = A*(D_S*G); FAD = FA*D_S;
        alpha_k = min(1,-(nabla_S(:)'*D_S(:))/...            
            (ADG(:)'*ADG(:)+FAD(:)'*FAD(:)+delta*D_S(:)'*D_S(:)));
        S = S_ex+alpha_k*D_S;
    else % (F)PG update
        nabla_S = (AtA*(S_ex*G)-AtY_H)*G'+FA'*(FA*S_ex-Y_M);
        if tight_flag == 1
            APhi = A*Phi;
        else
            APhi = A;
        end
        beta_S = max(norm(APhi'*thetaG_FFt*APhi,2),delta);
        S = Simplex_Proj(S_ex-1/beta_S*nabla_S);
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
    end
    
    % -------- check convergence 
    f_prev = f_curr;
    FA = F*A;
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
Y_HiBCD = A*S;
end

function maxEig = PowerMethod(M,iteraNum)
n = size(M,1); x = rand(n,1)-0.5; x = x./abs(x);
for k = 1 : iteraNum 
    x = M*x; 
end
maxEig = (x'*M*x)/(x'*x);
end

function S_Proj = Simplex_Proj(S)
S_Proj = max(bsxfun(@minus,S,max(bsxfun(@rdivide,cumsum(sort(S,1,'descend'),1)-1,(1:size(S,1))'),[],1)),0);
end