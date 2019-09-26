clc; clear; close all;
addpath([pwd,'/../methods/FUMI']);
addpath([pwd,'/../methods/CNMF']);
addpath([pwd,'/../methods/SupResPALM']);
addpath([pwd,'/../methods']);
addpath([pwd,'/../functions']);
addpath([pwd,'/data and specifications']);
% ============================================

load('Chikusei_1080.mat');
rng('default');
load('HS_spec.mat');
load('MS_spec.mat');
dataset = patch;
[Lx,Ly,M] = size(dataset);
L = Lx*Ly;

% --------generate data
maxIter = 3000;
stopThreshd = 1e-4;
N = 20;
dsRatio = 8;
GauSigma = 1.7;
kernelSize = 11;
[G,B] = Construct_Toeplitz_G(Lx,Ly,kernelSize,GauSigma,dsRatio);
F = Construct_F(HS_spec,MS_spec);

% -------set recorders
algNum = 5;
trial = 100;
SNRCase = 4;
runtime = zeros(algNum,trial,SNRCase);
iter_record = zeros(algNum,trial,SNRCase);
objRecord = cell(algNum,trial,SNRCase);
timeRecord = cell(algNum,trial,SNRCase);

for i = 1:trial
    % --------generate data
    Y = reshape(dataset,[],M)';
    Y_M = F*Y; Y_H = Y*G;
    V_H = randn(size(Y_H)); V_M = randn(size(Y_M));
    % --------different noisy scenarios
    for j = 2:1:2
        SNR = j*10;
        % ----observation generation
        YM_sigma = sqrt((sum(Y_M(:).^2)/(L))/(10^(SNR/10)))/sqrt(size(F,1));
        YM_noise = Y_M+YM_sigma*V_M;
        YH_sigma = sqrt((sum(Y_H(:).^2)/(L/dsRatio^2))/(10^(SNR/10)))/sqrt(size(F,2));
        YH_noise = Y_H+YH_sigma*V_H;
        HSI_noise = reshape(YH_noise',Lx/dsRatio,Ly/dsRatio,[]);
        MSI_noise = reshape(YM_noise',Lx,Ly,[]);
        
        % ---- initialization
        A_init = max(min(sisal(YH_noise,N, 'spherize', 'no','MM_ITERS',80, 'TAU', 0.006, 'verbose',0),1),0);
        S_tmp = sunsal(A_init,YH_noise,'POSITIVITY','yes','ADDONE','yes','verbose','no');
        S_init = Simplex_Proj(hyperConvert2d(imresize(hyperConvert3d(S_tmp,Lx/dsRatio),dsRatio)));
        obj_init = (norm(YH_noise-A_init*(S_init*G),'fro')^2+norm(YM_noise-F*A_init*S_init,'fro')^2)/2;
        fprintf('tiral: %g, SNR: %gdB:\n',i,SNR);
        Y_cut = cutBoundary(Y,Lx,Ly,kernelSize);
        
        % ---- HiBCD-FUMI (FW+FPG)
        caseCode = 1;
        timer = clock;
        [Y_FW,objRecord{caseCode,i,j},timeRecord{caseCode,i,j},iter_record(caseCode,i,j)] = HiBCD_stepSize(HSI_noise,MSI_noise,F,G,N,'S_UPDATE',2,'A_UPDATE',3,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init,'FW_STEP_SIZE',1);
        runtime(caseCode,i,j) = etime(clock,timer);
        
        % ---- HiBCD-FUMI (FW+FPG)
        caseCode = 2;
        timer = clock;
        [Y_FW,objRecord{caseCode,i,j},timeRecord{caseCode,i,j},iter_record(caseCode,i,j)] = HiBCD_stepSize(HSI_noise,MSI_noise,F,G,N,'S_UPDATE',2,'A_UPDATE',3,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init,'FW_STEP_SIZE',2);
        runtime(caseCode,i,j) = etime(clock,timer);
        
        % ---- HiBCD-FUMI (FW+FPG)
        caseCode = 3;
        timer = clock;
        [Y_FW,objRecord{caseCode,i,j},timeRecord{caseCode,i,j},iter_record(caseCode,i,j)] = HiBCD_stepSize(HSI_noise,MSI_noise,F,G,N,'S_UPDATE',2,'A_UPDATE',3,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init,'FW_STEP_SIZE',3);
        runtime(caseCode,i,j) = etime(clock,timer);
        
        % ---- HiBCD-FUMI (FW+FPG)
        caseCode = 4;
        timer = clock;
        [Y_FW,objRecord{caseCode,i,j},timeRecord{caseCode,i,j},iter_record(caseCode,i,j)] = HiBCD_stepSize(HSI_noise,MSI_noise,F,G,N,'S_UPDATE',3,'A_UPDATE',3,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init,'FPG_STEP_SIZE',1);
        runtime(caseCode,i,j) = etime(clock,timer);
        
        % ---- HiBCD-FUMI (FW+FPG)
        caseCode = 5;
        timer = clock;
        [Y_FW,objRecord{caseCode,i,j},timeRecord{caseCode,i,j},iter_record(caseCode,i,j)] = HiBCD_stepSize(HSI_noise,MSI_noise,F,G,N,'S_UPDATE',3,'A_UPDATE',3,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init,'FPG_STEP_SIZE',2);
        runtime(caseCode,i,j) = etime(clock,timer);
        
    end
    
    record.runtime = runtime;
    record.iter_record = iter_record;
    
end
