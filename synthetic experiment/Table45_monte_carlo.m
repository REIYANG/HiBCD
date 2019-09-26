clc; clear; close all;
addpath([pwd,'/../methods/FUMI']);
addpath([pwd,'/../methods/CNMF']);
addpath([pwd,'/../methods/SupResPALM']);
addpath([pwd,'/../methods']);
addpath([pwd,'/../functions']);
addpath([pwd,'/data and specifications']);
% ============================================
rng('default');
load('S.mat');
load('USGS_pruned.mat');
load('HS_spec.mat');
load('RGB_spec.mat');
load('MS_spec.mat');
Lx = 120; Ly = 120; L = Lx*Ly;
M = 224;
% --------generate data
maxIter = 3000;
stopThreshd = 1e-4;
N = 10;
dsRatio = 4;
GauSigma = 1.7;
kernelSize = 11;
[G,B] = Construct_Toeplitz_G(Lx,Ly,kernelSize,GauSigma,dsRatio);
F = Construct_F(HS_spec,MS_spec);
lx = Lx-kernelSize+1;
ly = Ly-kernelSize+1;

% ----CNMF setting
CNMF_opt.convthresh = stopThreshd;
CNMF_opt.maxIteraNum = maxIter;
CNMF_opt.F = F;
CNMF_opt.G = G;
% ----FUMI settings
FUMI_opt.M_M = size(F,1);
FUMI_opt.W1 = Lx;
FUMI_opt.W2 = Ly;
FUMI_opt.M = size(F,2);
FUMI_opt.dsRatio = dsRatio;
FUMI_opt.N = N;
FUMI_opt.F = F;
FUMI_opt.G = G;
FUMI_opt.GauSigma = GauSigma;
FUMI_opt.kernelSize = kernelSize;
FUMI_opt.maxIteration = maxIter;
FUMI_opt.convthresh = stopThreshd;
FUMI_opt.It_ADMM_EEA = 30; FUMI_opt.Th_ADMM_EEA = 1e-4;
FUMI_opt.It_ADMM_ABUN = 30; FUMI_opt.Th_ADMM_ABUN = 1e-4;
% ----HySure settings
HySure_opt.maxIteraNum = maxIter;
HySure_opt.F = F;
HySure_opt.G = G;
HySure_opt.Lx = Lx; HySure_opt.Ly = Ly;
HySure_opt.N = N;
HySure_opt.dsRatio = dsRatio;
HySure_opt.GauSigma = GauSigma;
HySure_opt.kernelSize = kernelSize;
HySure_opt.convthresh = stopThreshd;

% -------set recorders
algNum = 10;
trial = 100;
SNRCase = 4;
runtime = zeros(algNum,trial,SNRCase);
PSNR = zeros(algNum,trial,SNRCase);
SAM = zeros(algNum,trial,SNRCase);
ERGAS = zeros(algNum,trial,SNRCase);
iter_record = zeros(algNum,trial,SNRCase);
objRecord = cell(algNum,trial,SNRCase);
timeRecord = cell(algNum,trial,SNRCase);
for i = 1:trial
    % --------generate data
    Lx_pos = randi(size(ABUND,1)-Lx-1,1);
    Ly_pos = randi(size(ABUND,2)-Ly-1,1);
    S = reshape(ABUND(Lx_pos:Lx_pos+Lx-1,Ly_pos:Ly_pos+Ly-1,:),[],9)';
    A = USGS(:,randperm(size(USGS,2),9));
    Y = A*S;
    Y_M = F*Y; Y_H = Y*G;
    V_H = randn(size(Y_H)); V_M = randn(size(Y_M));
    % --------different noisy scenarios
    for j = 1:1:SNRCase
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
        
        % ---- naive interpretation
        caseCode = 1;
        timer = clock;
        img_resize = zeros(Lx,Ly,M);
        img_HS = reshape(YH_noise',Lx/dsRatio,Ly/dsRatio,[]);
        for band_count = 1:M
            img_resize(:,:,band_count) = imresize(img_HS(:,:,band_count),[Lx,Ly]);
        end
        runtime(caseCode,i,j) = etime(clock,timer);
        Y_interpolation = cutBoundary(reshape(img_resize,[],M)',Lx,Ly,kernelSize);
        [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)]  = ...
            quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_interpolation',lx,ly,[]),0,1/dsRatio);
        fprintf('Naive:  N=%g --- time: %g.\n',N,runtime(caseCode,i,j));
        fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
            PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
        
%         % ---- CNMF
%         caseCode = 2;
%         CNMF_opt.A_init = A_init;
%         CNMF_opt.S_init = S_init;
%         timer = clock;
%         [A_CNMF,S_CNMF,~,iter_record(caseCode,i,j)] = CNMF(HSI_noise,MSI_noise,CNMF_opt);
%         runtime(caseCode,i,j) = etime(clock,timer);
%         Y_CNMF = cutBoundary(A_CNMF*S_CNMF,Lx,Ly,kernelSize);
%         [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)] = ...
%             quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_CNMF',lx,ly,[]),0,1/dsRatio);
%         fprintf('CNMF: N=%g --- time: %g, iteration: %g.\n',N,runtime(caseCode,i,j),iter_record(caseCode,i,j));
%         fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
%             PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
%         
%         % ---- FUMI
%         caseCode = 3;
%         timer = clock;
%         FUMI_opt.A_init = A_init;
%         [FUMI_Out,para_opt]= FUMI(YH_noise,YM_noise,FUMI_opt);
%         runtime(caseCode,i,j) = etime(clock,timer);
%         timeRecord{caseCode,i,j} = para_opt.time_record;
%         Y_FUMI = cutBoundary(FUMI_Out.Y_FUMI,Lx,Ly,kernelSize);
%         iter_record(caseCode,i,j) = para_opt.iter;
%         [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)]  = ...
%             quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_FUMI',lx,ly,[]),0,1/dsRatio);
%         fprintf('FUMI: N=%g --- time: %g, iteration: %g.\n',N,runtime(caseCode,i,j),iter_record(caseCode,i,j));
%         fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
%             PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
%         
%         % ---- SupResPALM
%         caseCode = 4;
%         timer = clock;
%         [Y_SupResPALM,iter_record(caseCode,i,j),objRecord{caseCode,i,j},timeRecord{caseCode,i,j}] = SupResPALM(YH_noise,YM_noise,A_init,S_init,F,G,N,stopThreshd,maxIter);
%         runtime(caseCode,i,j) = etime(clock,timer);
%         Y_SupResPALM = cutBoundary(Y_SupResPALM,Lx,Ly,kernelSize);
%         [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)]  = ...
%             quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_SupResPALM',lx,ly,[]),0,1/dsRatio);
%         fprintf('SupResPALM: N=%g --- time: %g, iteration: %g.\n',N,runtime(caseCode,i,j),iter_record(caseCode,i,j));
%         fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
%             PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
        
        % ---- HiBCD-plain CoSMF (FW+FPG)
        caseCode = 5;
        timer = clock;
        [Y_FW,objRecord{caseCode,i,j},timeRecord{caseCode,i,j},iter_record(caseCode,i,j)] = HiBCD_plain(HSI_noise,MSI_noise,F,G,N,'S_UPDATE',2,'A_UPDATE',3,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init);
        runtime(caseCode,i,j) = etime(clock,timer);
        Y_plain_1 = cutBoundary(Y_FW,Lx,Ly,kernelSize);
        [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)]  = ...
            quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_plain_1',lx,ly,[]),0,1/dsRatio);
        fprintf('HiBCD-plain CoSMF(FW+FPG): N=%g --- time: %g, iteration: %g.\n',N,runtime(caseCode,i,j),iter_record(caseCode,i,j));
        fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
            PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
        %
        % ---- HiBCD-plain CoSMF (FPG+FPG)
        caseCode = 6;
        timer = clock;
        [Y_FW,objRecord{caseCode,i,j},timeRecord{caseCode,i,j},iter_record(caseCode,i,j)] = HiBCD_plain(HSI_noise,MSI_noise,F,G,N,'S_UPDATE',3,'A_UPDATE',3,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init);
        runtime(caseCode,i,j) = etime(clock,timer);
        Y_plain_2 = cutBoundary(Y_FW,Lx,Ly,kernelSize);
        [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)]  = ...
            quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_plain_2',lx,ly,[]),0,1/dsRatio);
        fprintf('HiBCD-plain CoSMF(FPG+FPG): N=%g --- time: %g, iteration: %g.\n',N,runtime(caseCode,i,j),iter_record(caseCode,i,j));
        fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
            PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
        
        % ---- HiBCD-plain CoSMF (FPG+FPG)
        caseCode = 7;
        timer = clock;
        [Y_FW,objRecord{caseCode,i,j},timeRecord{caseCode,i,j},iter_record(caseCode,i,j)] = HiBCD_plain(HSI_noise,MSI_noise,F,G,N,'S_UPDATE',2,'A_UPDATE',2,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init);
        runtime(caseCode,i,j) = etime(clock,timer);
        Y_plain_3 = cutBoundary(Y_FW,Lx,Ly,kernelSize);
        [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)]  = ...
            quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_plain_3',lx,ly,[]),0,1/dsRatio);
        fprintf('HiBCD-plain CoSMF(FW+FW): N=%g --- time: %g, iteration: %g.\n',N,runtime(caseCode,i,j),iter_record(caseCode,i,j));
        fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
            PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
        
        % ---- HiBCD-NNC CoSMF (FW+FPG)
        caseCode = 8;
        timer = clock;
        [Y_PG,objRecord{caseCode,i,j,:},timeRecord{caseCode,i,j,:},iter_record(caseCode,i,j)] = HiBCD_NNC(HSI_noise,MSI_noise,F,G,N,'A_UPDATE',3,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init,'Gamma',10,'S_UPDATE',2);
        runtime(caseCode,i,j) = etime(clock,timer);
        Y_nuclear_1 = cutBoundary(Y_PG,Lx,Ly,kernelSize);
        [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)]  = ...
            quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_nuclear_1',lx,ly,[]),0,1/dsRatio);
        fprintf('HiBCD-NNC CoSMF(FW+FPG): N=%g --- time: %g, iteration: %g.\n',N,runtime(caseCode,i,j),iter_record(caseCode,i,j));
        fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
            PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
        
        % ---- HiBCD-NNC CoSMF (FPG+FPG)
        caseCode = 9;
        timer = clock;
        [Y_PG,objRecord{caseCode,i,j,:},timeRecord{caseCode,i,j,:},iter_record(caseCode,i,j)] = HiBCD_NNC(HSI_noise,MSI_noise,F,G,N,'A_UPDATE',3,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init,'Gamma',10,'S_UPDATE',3);
        runtime(caseCode,i,j) = etime(clock,timer);
        Y_nuclear_2 = cutBoundary(Y_PG,Lx,Ly,kernelSize);
        [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)]  = ...
            quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_nuclear_2',lx,ly,[]),0,1/dsRatio);
        fprintf('HiBCD-NNC CoSMF(FPG+FPG): N=%g --- time: %g, iteration: %g.\n',N,runtime(caseCode,i,j),iter_record(caseCode,i,j));
        fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
            PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
        
        % ---- HiBCD-NNC CoSMF (FW+FW)
        caseCode = 10;
        timer = clock;
        [Y_PG,objRecord{caseCode,i,j,:},timeRecord{caseCode,i,j,:},iter_record(caseCode,i,j)] = HiBCD_NNC(HSI_noise,MSI_noise,F,G,N,'A_UPDATE',2,...
            'MAXITER',maxIter,'STOPPING',stopThreshd,'S_INIT',S_init,'A_INIT',A_init,'Gamma',10,'S_UPDATE',2);
        runtime(caseCode,i,j) = etime(clock,timer);
        Y_nuclear_3 = cutBoundary(Y_PG,Lx,Ly,kernelSize);
        [PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j)]  = ...
            quality_assessment(reshape(Y_cut',lx,ly,[]),reshape(Y_nuclear_3',lx,ly,[]),0,1/dsRatio);
        fprintf('HiBCD-NNC CoSMF(FW+FW): N=%g --- time: %g, iteration: %g.\n',N,runtime(caseCode,i,j),iter_record(caseCode,i,j));
        fprintf('               PSNR: %g, ERGAS: %g, SAM: %g.\n',...
            PSNR(caseCode,i,j),ERGAS(caseCode,i,j),SAM(caseCode,i,j));
    end
    record.runtime = runtime;
    record.iter_record = iter_record;
    record.timeRecord = timeRecord;
    record.objRecord = objRecord;
    record.PSNR = PSNR;
    record.SAM = SAM;
    record.ERGAS = ERGAS;
end
