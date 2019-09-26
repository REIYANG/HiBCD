% Code (c) Charis Lanaras, ETH Zurich, Oct 28 2015 (see LICENSE)
% charis.lanaras@geod.baug.ethz.ch

function [Y,j,f_record,time_record] = SupResPALM(hyper, multi, E_init,A_init,F, G, p,convThresh,maxIter )
% SupResPALM - Perform hyperspectral super-resolution by spectral unmixing

% Usage
%   [E,A] = SupResPALM(hyper, multi, truth, srf, p, h )
% 
% Inputs
%   hyper - hyperspectral image (in 2D format)
%   multi - RGB or multispectral image (in 2D format)
%   truth - ground truth image (in 2D format) - used only for evaluation
%   srf   - the spectral response function of the RGB multispectral camera
%   p     - the desired number of material spectra to extract
%   h     - (optional) the hight of the image, in case of non-square images
%           (h must be in accordance with the spatial downsampling factor)
%
% Outputs
%   E - Matrix of endmembers (spectral basis)
%   Y - Abundances (mixing coefficients)
%
% References
%   C. Lanaras, E. Baltsavias, K. Schindler. Hyperspectral Super-Resolution
%   by Coupled Spectral Unmixing. In: ICCV 2015
% 
% Comment
%    To transform a 3D image cube to the respective 2D format, you can use
%    the hyperConvert2d.m function.
%

if ndims(hyper)==3
    hyper = hyperConvert2d(hyper);
end
if ndims(multi)==3
    multi = hyperConvert2d(multi);
end
% if ndims(truth)==3
%     truth = hyperConvert2d(truth);
% end
if ~(nargin==6)
    h = sqrt(size(multi,2));
end

epsilon = 0.0001;
s2=p/2; %Active endmemebers per pixel (on average over the whole image)
% Default value
% To activate this option you need to uncomment the specified line in
%   highResStep.m
resAB(1) = 1;

scale_diff = sqrt(size(multi,2)/size(hyper,2)); % difference of resolution
% [S St] = hyperSpatialDown(h, size(multi,2)/h, scale_diff);

% Initialisations
E = E_init;
% sisal(hyper,p, 'spherize', 'no','MM_ITERS',80, 'TAU', 0.006, 'verbose',0);
% AS = sunsal(E, hyper,'POSITIVITY','yes','ADDONE','yes','verbose','no');
A = A_init;
AS = A*G;
% hyperConvert2d(imresize(hyperConvert3d(AS,h/scale_diff),scale_diff));

psnr_record = [];
sam_record = [];
obj_curr = 0;
time_record = zeros(maxIter+1,1);
f_record = zeros(maxIter+1,1);
tic;
f_curr = 0.5 * sum(sum((hyper-E*A*G).^2)) + 0.5 * sum(sum((multi-F*E*A).^2));
f_record(1) = f_curr;
for j=1:maxIter
    
    % hyperspectral least-squares
    [E, ~, res] = lowResStep(hyper,E,AS);
    resA2(j) = min(res);
    
    % spectraly downgrade endmemebrs
    RE = F*E;
    
    % multispectral least-squares
    [~, A, res] = highResStep(multi,RE,A,s2);
    resB2(j) = min(res);
    
    % update abundances for hyperspectral step
    AS = A*G;
    
    % Residual of the objective function (5a)
    resAB(j+1) = resA2(j)+resB2(j);
    
    % Compute RMSE only for printing during procedure
%     RMSE(j) = hyperErrRMSE(truth,E*A);
    
%     % Convergence checks
%     if ( resAB(j) / resAB(j+1) ) > 1+epsilon || ( resAB(j) / resAB(j+1) ) < 1-epsilon
%         fprintf(['Iter: ' num2str(j) ' RMSE: ' num2str(RMSE(j)) '\n'])
%     else
%         fprintf(['Iter: ' num2str(j) ' RMSE: ' num2str(RMSE(j)) '\n'])
%         fprintf(['Stopped after ' num2str(j) ' iterations. Final RMSE: ' num2str(RMSE(j)) '\n'])
%         break
%     end
%     [psnr_tmp,sam_tmp] = quality_mini(reshape(Y',W1,W2,[]),reshape((E*A)',W1,W2,[]),0);
%         psnr_record = [psnr_record psnr_tmp];
%         sam_record = [sam_record sam_tmp];
%         subplot(2,1,1); plot(psnr_record,'LineWidth',3);
%         subplot(2,1,2); plot(sam_record,'LineWidth',3);pause(0.001)
    
    obj_prev = obj_curr ;
    obj_curr = 0.5 * sum(sum((hyper-E*A*G).^2)) + 0.5 * sum(sum((multi-F*E*A).^2)) ;
    time_record(j+1) = toc+time_record(j);
    tic;
    f_record(j+1) = obj_curr;
    obj_chng = abs(obj_curr-obj_prev)/obj_curr ;
    if( obj_chng < convThresh )
        break ;
    end
end
f_record(j+2:end) = [];
time_record(j+2:end) = [];
Y = E*A;
end

function [ E, A, res ] = highResStep( M, E, A, sparse_factor )
% Solving eq. (7) with a projected gradient descent method.

maxIter = 100;
epsilon = 1.01;
gamma2 = 1.01;

N = size(M, 2);
beta = round(sparse_factor*N); % Number of desirable non-zero entries

res(1) = norm(M-E*A,'fro')+100;

for k=1:maxIter
    E_old = E;
    A_old = A;
    
    % 2.2. Update the Abundances
    dk = gamma2 * norm( E*E' ,'fro');
    V = A - 1/dk * E' * ( E*A - M );
    
    % Uncomment Tau_multi and comment the following line to use the sparse
    % constraint
    % A = Tau_multi(Pplusb(V),beta);
    A = Pplusb(V);
    
    % Calculation of residuals
    res(k+1,1) = sqrt(norm(M-E*A,'fro')^2/size(M,1)/size(M,2));
    
    % Checks for exiting iteration
    if (1/res(k+1) * res(k)) < epsilon
%         fprintf(['Multi: Iter ' num2str(k) ', res: ' num2str(res(k+1)*1000) '. '])
        break
    end
    
    if (res(k+1) / res(k))>1
        E = E_old;
        A = A_old;
%         text = ['Multi: Exited during ' num2str(k) 'th iteration with residual ' num2str(res(k+1)*1000) ', because res increased'];
%         disp(text)
        break
    end
    
end
end

function [ E, A, res ] = lowResStep( H, E, A )
% Solving eq. (6) with a projected gradient descent method.

maxIter = 100;
epsilon = 1.01;
gamma1 = 1.01;

notfirst = 0;

res(1) = norm(H-E*A,'fro')+100;

for k=1:maxIter
    E_old = E;
    A_old = A;
    
    % 2.1. Update of signatures
    ck = gamma1 * norm(A*A','fro');
    U = E - 1/ck * ( E*A - H ) * A';
    E = Pplusa(U);
    
    % Calculation of residuals
    res(k+1,1) = sqrt(norm(H-E*A,'fro')^2/size(H,1)/size(H,2));
    
    % Checks for exiting iteration
    if (1/res(k+1) * res(k)) < epsilon
%         fprintf(['Hyper: Iter ' num2str(k) ', res: ' num2str(res(k+1)*1000) '. '])
        break
    end
    
    if (res(k+1) / res(k))>1
        if notfirst == 1
            E = E_old;
            A = A_old;
%             text = ['Hyper: Exited during ' num2str(k) 'th iteration with residual ' num2str(res(k+1)*1000) ', because res increased'];
            disp(text)
            break
        else
            notfirst = 1;
        end
    end
    
end
end

function U = Pplusa(U)
% max{0,U}
U(U<0) = 0;
U(U>1) = 1;
end

function V = Pplusb(V)
% Simplex Projection

V = Simplex_Proj(V);

end

function U = Tau_multi(U,s)

% keep only the first s largest entries of U
U1 = reshape(U,[],1);
[values, ind] = sort(U1,'descend');
U1 = zeros(length(U1),1);
U1(ind(1:s),1) = values(1:s);
U = reshape(U1,size(U));

end
