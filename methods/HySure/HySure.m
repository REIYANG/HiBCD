function [ Zimhat,f_record,time_record] = HySure( Yhim, Ymim, downsamp_factor, R, B, p, basis_type, lambda_phi, lambda_m, Setting )
%data_fusion - The actual data fusion algorithm
% 
%   Given the two observed images (HS + PAN/MS), fuses the data according
%   to the steps described in detail in the Appendix of [1]. 
% 
% [ Zimhat ] = data_fusion( Yhim, Ymim, downsamp_factor, R, B, p, 
%                             basis_type, lambda_phi, lambda_m )
% 
% Input: 
% Yhim: HS image, 
% Ymim: MS image, 
% downsamp_factor: downsampling factor,
% R: relative spectral resppnse;
% B: relative spatial response;
% p: corresponds to variable L_s in [1]; number of endmembers in VCA /
%     number of non-truncated singular vectors,
% basis_type: method to estimate the subspace. Can be either 'VCA' or
%     'SVD'
% lambda_phi: regularization parameter lambda_phi (corresponds to VTV), 
% lambda_m: regularization parameter lambda_m.
% 
% Output: 
% Zimhat: estimated image with high spatial and spectral resolution
% 
%   [1] M. Simoes, J. Bioucas-Dias, L. Almeida, and J. Chanussot, 
%        ?A convex formulation for hyperspectral image superresolution 
%        via subspace-based regularization,? IEEE Trans. Geosci. Remote 
%        Sens., to be publised.

% % % % % % % % % % % % % 
% 
% Version: 1
% 
% Can be obtained online from: https://github.com/alfaiate/HySure
% 
% % % % % % % % % % % % % 
% 
% Copyright (C) 2014 Miguel Simoes, Jose Bioucas-Dias, Luis B. Almeida 
% and Jocelyn Chanussot
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% 
% % % % % % % % % % % % % 
% 
% This funtion has three steps:
% 
% I. Precomputations (for example, FFTs of the filters that will be used 
% throughout the function). 
% 
% II. It will then learn the subspace where the HS ''lives'', via SVD
% or VCA.
% 
% Basis type - 1: find endmembers with VCA, 3: learn the
%     subspace with SVD
% 
% III. The iterative process used to compute a solution to the
% optimization problem via ADMM/SALSA. The following parameters can be 
% adjusted:
% ADMM parameter
mu = 0.05;
% ADMM iterations
iters = Setting.maxIteraNum;
G = Setting.G;
time_record = zeros(iters+1,1);
f_record = zeros(iters+1,1);
tic
% 1
% IV. Postprocessing (denoising).
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% I. Precomputations.                                                   %
% -------------------                                                   %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

[nl, nc, ~] = size(Ymim);
% Define the difference operators as filters
dh = zeros(nl,nc);
dh(1,1) = 1;
dh(1,nc) = -1;

dv = zeros(nl,nc);
dv(1,1) = 1;
dv(nl,1) = -1;

FDH = fft2(dh);
FDHC = conj(FDH);
FDV = fft2(dv);
FDVC = conj(FDV);

% Fourier transform of B

B = double(B);
FB = fft2(B);
FBC = conj(FB);

IBD_B  = FBC  ./(abs(FB.^2) + abs(FDH).^2+ abs(FDV).^2 + 1);
IBD_II = 1    ./(abs(FB.^2) + abs(FDH).^2+ abs(FDV).^2 + 1);
IBD_DH = FDHC ./(abs(FB.^2) + abs(FDH).^2+ abs(FDV).^2 + 1);
IBD_DV = FDVC ./(abs(FB.^2) + abs(FDH).^2+ abs(FDV).^2 + 1);

% % % % % % % % % % % %
% We will work with image Yhim in a matrix with the same size as Ymim. The
% result is a matrix filled with zeros. We do this for computational 
% convenience and the end result is the same. We also explicity form a
% subsampling mask that has the same effect has matrix M in [1].
shift = 1;
mask = zeros(nl, nc);
mask(shift+1:downsamp_factor:nl, shift+1:downsamp_factor:nc) = 1;
% Subsampling mask in image format
maskim = repmat(mask, [1, 1, p]);
% Subsampling mask in matrix format
mask = im2mat(maskim);
% Yhim with the same size as Ym (but with zeros)
Yhim_up = upsamp_HS(Yhim, downsamp_factor, nl, nc, shift);
Yh_up = im2mat(Yhim_up);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% II. Subspace learning.                                                %
% ----------------------                                                %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
switch  basis_type 
    case 'VCA'
%     Find endmembers with VCA (pick the one with smallest volume from 20 
%     runs of the algorithm)
    max_vol = 0;
    vol = zeros(1, 20);
    for idx_VCA = 1:20
        E_aux = VCA(Yh_up(:, mask(1,:) > 0),'Endmembers',p,'SNR',0,'verbose','off');
        vol(idx_VCA) = abs(det(E_aux'*E_aux));
        if vol(idx_VCA) > max_vol
            E = E_aux;
            max_vol = vol(idx_VCA);
        end   
    end
    
    case   basis_type == 'SVD'
%     Learn the subspace with SVD
    Ry = Yh(:, mask(1,:) > 0)*Yh(:, mask(1,:) > 0)'/np;
    Ry = Yh;
    [E, ~] = svds(Ry,p);

end
A = E;
[W1,W2,~] = size(Ymim);
FA = R*A;
D_v = kron(sparse(eye(W2)),sparse(gallery('circul',[1 -1 zeros(1,W1-2)])));
D_h = kron(sparse(gallery('circul',[1 -1 zeros(1,W2-2)])),sparse(eye(W1)));
rho = lambda_phi;
Y_H = reshape(Yhim,[],size(Yhim,3))';
Y_M = reshape(Ymim,[],size(Ymim,3))';
convthresh = Setting.convthresh;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% III. ADMM/SALSA.                                                      %
% ----------------                                                      %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Auxiliary matrices
p = size(E,2);
mask(p+1:end,:) = [];
IE = E'*E+mu*eye(p);
yyh = E'*Yh_up;

IRE = lambda_m*E'*(R'*R)*E+mu*eye(p);
Ym = im2mat(Ymim);
yym = E'*R'*Ym;

% Define and initialize variables
% All variables (X,V1,V2,V3,D1,D2,D3) will be represented as matrices
X = zeros(nl*nc,p)';
V1 = X;
D1 = X;
V2 = X;
D2 = X;
V3 = X;
D3 = X;
V4 = X;
D4 = X;
objPrev = 10;
% costf = [];

% Initialize tic
% t1 = tic; t2 = toc(t1); t1 = tic; t2 = toc(t1);
% t1 = tic;

f_record(1) = obj_Tol(A,FA,X,G,Y_H,Y_M,D_h,D_v,rho);
for i=1:iters
    
    %   min   ||XB - V1 - A1||_F^2  +
    %    X    ||X  - V2 - A2||_F^2
    %         ||XDH - V3 - A3||_F^2 +
    %         ||XDV - V4 - A4||_F^2
    %
    X = ConvC(V1+D1, IBD_B, nl) + ConvC(V2+D2, IBD_II, nl) + ...
        ConvC(V3+D3, IBD_DH, nl) +  ConvC(V4+D4, IBD_DV, nl);
    
    
    %  max (1/2)||Yh - EV1M|_F^2 + (mu/2)||XB - V1 - D1||_F^2
    %   V1
    NU1 =  ConvC(X, FB, nl) - D1;
    V1 = IE\(yyh + mu*NU1).*mask + NU1.*(1-mask);
    
    
    %  max (lambda_m/2)||Ym - REV2|_F^2 + (mu/2)||X - V2 - D2||_F^2
    %   V1
    NU2 =  X - D2;
    V2 = IRE\(lambda_m*yym + mu*NU2);
    
    
    % min lambda VTV(V2,V3) + (mu/2)||XDH - V3 - D3||_F^2 + (mu/2)||XDV - V4 - D4||_F^2
    % V2,V3
    NU3 =  ConvC(X, FDH, nl) - D3;
    NU4 =  ConvC(X, FDV, nl) - D4;
    [V3,V4] = vector_soft_col_iso(NU3,NU4,lambda_phi/mu);
    
%     fprintf('iter = %d out of %d\n', i, iters);
        
    % Update Lagrange multipliers
    D1 = -NU1 + V1;    % D1 - (XB - V1)
    D2 = -NU2 + V2;    % D2 - (XDH - V2)
    D3 = -NU3 + V3;    % D3 - (XDV - V3)
    D4 = -NU4 + V4;    % D3 - (XDV - V3)
    % ======================== %
    %   stopping criterion 1   %
    % ======================== %
    objCurr = obj_Tol(A,FA,X,G,Y_H,Y_M,D_h,D_v,rho);
    time_record(i+1) = toc+time_record(i);
    tic;
    f_record(i+1) = objCurr;
    obj_chng = abs(objPrev-objCurr)/objPrev;
    objPrev = objCurr;
    if obj_chng<convthresh && i>50
        break ;
    end
end
f_record(i+2:end) = [];
time_record(i+2:end) = [];
% t2 = toc(t1);
%fprintf('The algorithm took %2.2f seconds to run.\n', t2);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% IV. Postprocessing.                                                   %
% -------------------------------                                       %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

Zhat = E*X;
Zimhat = mat2im(Zhat, nl);
end

function [objS] = obj_Tol(A,FA,S,G,Y_H,Y_M,D_h,D_v,rho)
    objS = 1/2*sum(sum((Y_M-FA*S).^2))+1/2*sum(sum((Y_H-A*(S*G)).^2))+rho*sum((sum((S*D_h).^2)+sum((S*D_v).^2)+eps).^(1/2));
end

function [Y1,Y2] = vector_soft_col_iso(X1,X2,tau)
%
%  computes the isotropic vector soft columnwise


NU = sqrt(sum(X1.^2)+sum(X2.^2));
A = max(0, NU-tau);
A = repmat((A./(A+tau)),size(X1,1),1);
Y1 = A.*X1;
Y2 = A.*X2;
end

function [Ae, indice, Rp] = VCA(R,varargin)

% Vertex Component Analysis Algorithm [VCA]
%
% [VCA] J. Nascimento and J. Bioucas-Dias
% "Vertex component analysis: a fast algorithm to unmix hyperspectral data"
% IEEE Transactions on Geoscience and Remote Sensing,  vol. 43, no. 4, 
% pp. 898-910, 2005.
%
% -------------------------------------------------------------------
% Usage:
%
% [Ae, indice, Rp ]= vca(R,'Endmembers',p,'SNR',r,'verbose',v)
%
% ------- Input variables -------------------------------------------
%
%  R - matrix with dimensions L(channels) x N(pixels)
%      Each pixel is a linear mixture of p endmembers
%      signatures R = M X, where M  and X are the mixing matrix 
%      and the abundance fractions matrix, respectively.
%
% 'Endmembers'
%          p - number of endmembers in the scene
%
% ------- Output variables -------------------------------------------
%
% A      - estimated mixing matrix (endmembers signatures)
%
% indice - pixels chosen to be the most pure
%
% Rp     - Data R projected on the identified signal subspace
%
% ------- Optional parameters -----------------------------------------
%
% 'SNR'     - (double) signal to noise ratio (dB)
%             SNR is used to decide the type of projection: projective
%             or orthogonal.
%
% 'verbose' - [{'on'} | 'off']
% ---------------------------------------------------------------------
%
% Please see [VCA] for more details or contact the Authors
%
% -----------------------------------------------------------------------
% version: 3.0 (21-January-2012)
%
% Modifications w.r.t. version 2.1:
%     
%  - Increased efficiency in the memory usage
%  - Correction of a bug in SNR estimation
%  - detection of outliers in the projective projection
%
% -----------------------------------------------------------------------
% Copyright (2012):  Jos??? Nascimento (zen@isel.pt)
%                    Jos??? Bioucas Dias (bioucas@lx.it.pt)
%
% affineProj is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = 'on'; % default
snr_input = 0;  % estimate the SNR
p = 0;          % default number of endmembers


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim_in_par = length(varargin);
if (nargin - dim_in_par)~=1
    error('Wrong parameters');
elseif rem(dim_in_par,2) == 1
    error('Optional parameters should always go by pairs');
else
    for i = 1 : 2 : (dim_in_par-1)
        switch lower(varargin{i})
            case 'verbose'
                verbose = varargin{i+1};
            case 'endmembers'
                p = varargin{i+1};
            case 'snr'
                SNR = varargin{i+1};
                snr_input = 1;       % user input SNR
            otherwise
                fprintf(1,'Unrecognized parameter:%s\n', varargin{i});
        end %switch
    end %for
end %if


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(R)
    error('there is no data');
else
    [L N]=size(R);  % L number of bands (channels)
    % N number of pixels (LxC)
end

if (p<=0 | p>L | rem(p,1)~=0),
    error('ENDMEMBER parameter must be an  integer between 1 and L');
end

if (L-p < p) & (snr_input == 0)
    if strcmp (verbose, 'on'),
        fprintf(1,' i can not  estimate SNR [(no bands)-p < p]\n');
        fprintf(1,' i will apply the projective projection\n');
        snr_input = 1;
        SNR = 100;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


r_m = mean(R,2);
Corr = R*R'/N;
[Ud,Sd] = svds(Corr,p);            % computes the a p-orth basis


if snr_input == 0
    % estimate SNR
    Pt = trace(Corr);      % total power     (signal + noise)
    Pp = sum(diag(Sd));   % projected power (signal + projected noise)
    Pn = (Pt-Pp)/(1-p/L);       % rough noise power estimate
    if Pn > 0
        SNR = 10*log10(Pp/Pn);  % signal-to-noise ratio in dB
        if strcmp (verbose, 'on'), 
            fprintf(1,'SNR estimated = %g[dB]\n',SNR); 
        end
    else
        SNR = 0;
        if strcmp (verbose, 'on'), 
            fprintf(1,'Input data belongs to a p-subspace'); 
        end
    end
end

% SNR threshold to decide the projection:
%       Projective Projection
%       Projection on the  p-1 subspace
SNR_th = 15 + 10*log10(p);

if SNR < SNR_th,
    if strcmp (verbose, 'on')
        fprintf(1,'Select proj. on the to (p-1) subspace.\n')
        fprintf(1,'I will apply the projective projection\n')
    end
    % orthogonal projection on an the best (p-1) affine set 
    d = p-1;
    Cov  = Corr - r_m*r_m';
    [Ud,Sd] = svds(Cov,d);         % computes the a d-orth basis
    R_o = R - repmat(r_m,[1 N]);   % remove mean 
    x_p =  Ud' * R_o;  % project the zero-mean data onto  a p-subspace
    
    %  original data projected on the indentified subspace 
    if (d > size(x_p,1)) 
        Rp =  Ud * x_p(:,:) + repmat(r_m,[1 N]);   % again in dimension L
    else
        Rp =  Ud * x_p(1:d,:) + repmat(r_m,[1 N]);   % again in dimension L
    end
    % compute the angles (cosines) between the projected vectors and the
    % original
    cos_angles = sum(Rp.*R)./(sqrt(sum(Rp.^2).*sum(R.^2)));
    
    
    % lift to p-dim
    c = max(sum(x_p.^2,1))^0.5;
    y = [x_p ; c*ones(1,N)] ;
    
else
    if strcmp (verbose, 'on'), 
        fprintf(1,'... Select the projective proj. (dpft)\n');
    end
    
    % projective projection
    d = p;    
    % project into a p-dim subspace (filter noise)
    x_p = Ud'*R;

    %  original data projected on the indentified subspace 
    Rp =  Ud * x_p(1:d,:);      % again in dimension L 
    
    
    
    % find a direction orthogonal to the affine set
    u = mean(x_p,2)*p;             
    
    % ensure that angle(xp(:,i),u) is positive

    scale = sum( x_p .* repmat(u,[1 N]) );

    th = 0.01;   % close to zero
    mask = scale < th ; 
    scale = scale.*(1-mask) + mask;
    
    y =  x_p./ repmat(scale,[d 1]) ;
    pt_errors = find(mask);
    
    % replace the bad vectors with a vector in the middle of the simplex
    y(:,pt_errors) = (u/norm(u)^2)*ones(1,length(pt_errors));
    
     
end

%%%%%%%%%%%%%%%%%%%%%%%
% VCA algorithm
%%%%%%%%%%%%%%%%%%%%%%%
p = size(y,1);
% pointers for endmembers
indice = zeros(1,p);
% save endmembers
A = zeros(p,p);
A(p,1) = 1;

for i=1:p
    w = rand(p,1);
    f = w - A*pinv(A)*w;
    f = f / sqrt(sum(f.^2));
    v = f'*y;
    [v_max indice(i)] = max(abs(v));
    A(:,i) = y(:,indice(i));        % same as x(:,indice(i))
end
Ae = Rp(:,indice);

return;
%%%%%%%%%%%%%%%%%%%%%%%
% End of the vca function
%%%%%%%%%%%%%%%%%%%%%%%



end

function [Yhim_up] = upsamp_HS(Yhim, downsamp_factor, nl, nc, shift)
%upsamp_HS - convert image Ymim to an image matrix with the same size as 
% Ymim. The result is a matrix filled with zeros.

[nlh, nch, L] = size(Yhim);
aux = zeros(nlh*downsamp_factor, nch*downsamp_factor, L);
for i=1:L
    aux(:,:,i) = upsample(upsample(Yhim(:,:,i), downsamp_factor, shift)', downsamp_factor, shift)';
end
Yhim_up = aux(1:nl, 1:nc, :);
end

function [A] = mat2im(X, nl)
%mat2im - converts a matrix to a 3D image

[p, n] = size(X);
nc = n/nl;
A = reshape(X', nl, nc, p);
end

function [A] = im2mat(X)
%im2mat - converts a 3D image to a matrix

[nl, nc, p] = size(X);
A = reshape(X, nl*nc, p)';
end

function [Yhim] = downsamp_HS(Yhim_up, downsamp_factor, shift)
%downsamp_HS - the equivalent of applying matrix M

[nl, nc, L] = size(Yhim_up);
Yhim = zeros(ceil(nl/downsamp_factor), ceil(nc/downsamp_factor), L);
for i=1:L
    Yhim(:,:,i) = downsample(downsample(Yhim_up(:,:,i), downsamp_factor, shift)', downsamp_factor, shift)';
end
end

function [A] = ConvC(X, FK, nl)
%ConvC - defines a circular convolution (the same for all bands) accepting 
% a matrix and returnig a matrix. FK is the fft2 of a one-band filter 

[p, n] = size(X);
nc = n/nl;
A = reshape(real(ifft2(fft2(reshape(X', nl, nc, p)).*repmat(FK,[1, 1, p]))), nl*nc, p)';
end

function [quality, quality_map] = img_qi(img1, img2, block_size)

%========================================================================
%
%Copyright (c) 2001 The University of Texas at Austin
%All Rights Reserved.
% 
%This program is free software; you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation; either version 2 of the License, or
%(at your option) any later version.
% 
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
% 
%The GNU Public License is available in the file LICENSE, or you
%can write to the Free Software Foundation, Inc., 59 Temple Place -
%Suite 330, Boston, MA 02111-1307, USA, or you can find it on the
%World Wide Web at http://www.fsf.org.
%
%Author  : Zhou Wang 
%Version : 1.0
% 
%The authors are with the Laboratory for Image and Video Engineering
%(LIVE), Department of Electrical and Computer Engineering, The
%University of Texas at Austin, Austin, TX.
%
%Kindly report any suggestions or corrections to zwang@ece.utexas.edu
%
%Acknowledgement:
%The author would like to thank Mr. Umesh Rajashekar, the Matlab master
%in our lab, for spending his precious time and giving his kind help
%on writing this program. Without his help, this program would not
%achieve its current efficiency.
%
%========================================================================
%
%This is an efficient implementation of the algorithm for calculating
%the universal image quality index proposed by Zhou Wang and Alan C. 
%Bovik. Please refer to the paper "A Universal Image Quality Index"
%by Zhou Wang and Alan C. Bovik, published in IEEE Signal Processing
%Letters, 2001. In order to run this function, you must have Matlab's
%Image Processing Toobox.
%
%Input : an original image and a test image of the same size
%Output: (1) an overall quality index of the test image, with a value
%            range of [-1, 1].
%        (2) a quality map of the test image. The map has a smaller
%            size than the input images. The actual size is
%            img_size - BLOCK_SIZE + 1.
%
%Usage:
%
%1. Load the original and the test images into two matrices
%   (say img1 and img2)
%
%2. Run this function in one of the two ways:
%
%   % Choice 1 (suggested):
%   [qi qi_map] = img_qi(img1, img2);
%
%   % Choice 2:
%   [qi qi_map] = img_qi(img1, img2, BLOCK_SIZE);
%
%   The default BLOCK_SIZE is 8 (Choice 1). Otherwise, you can specify
%   it by yourself (Choice 2).
%
%3. See the results:
%
%   qi                    %Gives the over quality index.
%   imshow((qi_map+1)/2)  %Shows the quality map as an image.
%
%========================================================================

if (nargin == 1 | nargin > 3)
   quality = -Inf;
   quality_map = -1*ones(size(img1));
   return;
end

if (size(img1) ~= size(img2))
   quality = -Inf;
   quality_map = -1*ones(size(img1));
   return;
end

if (nargin == 2)
   block_size = 8;
end

N = block_size.^2;
sum2_filter = ones(block_size);

img1_sq   = img1.*img1;
img2_sq   = img2.*img2;
img12 = img1.*img2;

img1_sum   = filter2(sum2_filter, img1, 'valid');
img2_sum   = filter2(sum2_filter, img2, 'valid');
img1_sq_sum = filter2(sum2_filter, img1_sq, 'valid');
img2_sq_sum = filter2(sum2_filter, img2_sq, 'valid');
img12_sum = filter2(sum2_filter, img12, 'valid');

img12_sum_mul = img1_sum.*img2_sum;
img12_sq_sum_mul = img1_sum.*img1_sum + img2_sum.*img2_sum;
numerator = 4*(N*img12_sum - img12_sum_mul).*img12_sum_mul;
denominator1 = N*(img1_sq_sum + img2_sq_sum) - img12_sq_sum_mul;
denominator = denominator1.*img12_sq_sum_mul;

quality_map = ones(size(denominator));
index = (denominator1 == 0) & (img12_sq_sum_mul ~= 0);
quality_map(index) = 2*img12_sum_mul(index)./img12_sq_sum_mul(index);
index = (denominator ~= 0);
quality_map(index) = numerator(index)./denominator(index);

quality = mean2(quality_map);
end