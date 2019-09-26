clc; clear; close all;
addpath([pwd,'/data and specifications']);
load('record.mat');
load('HS_spec.mat');
colorscheme = distinguishable_colors(9);
PSNR_curve = record.psnr_curve;
% ---- add markers
lw = 3;
ms = 20;
mint = 30;
plot(HS_spec(1:mint:end),PSNR_curve(1,1:mint:end),'s','MarkerSize',ms,'Linewidth',lw,'Color',colorscheme(1,:));
hold on;
plot(HS_spec(1:mint:end),PSNR_curve(2,1:mint:end),'o','MarkerSize',ms,'Linewidth',lw,'Color',colorscheme(7,:));
plot(HS_spec(1:mint:end),PSNR_curve(3,1:mint:end),'+','MarkerSize',ms,'Linewidth',lw,'Color',colorscheme(8,:));
plot(HS_spec(1:mint:end),PSNR_curve(4,1:mint:end),'p','MarkerSize',ms,'Linewidth',lw,'Color',colorscheme(4,:));
plot(HS_spec(1:mint:end),PSNR_curve(6,1:mint:end),'>','MarkerSize',ms,'Linewidth',lw,'Color',colorscheme(9,:));
plot(HS_spec(1:mint:end),PSNR_curve(12,1:mint:end),'v','MarkerSize',ms,'Linewidth',lw,'Color',colorscheme(6,:));
% ---- plot complete curves
int = 2;
plot(HS_spec(1:int:end),PSNR_curve(1,1:int:end),'-.','Linewidth',lw,'Color',colorscheme(1,:));
hold on;
plot(HS_spec(1:int:end),PSNR_curve(2,1:int:end),'-.','Linewidth',lw,'Color',colorscheme(7,:));
plot(HS_spec(1:int:end),PSNR_curve(3,1:int:end),'-.','Linewidth',lw,'Color',colorscheme(8,:));
plot(HS_spec(1:int:end),PSNR_curve(4,1:int:end),'-.','LineWidth',lw,'Color',colorscheme(4,:));
plot(HS_spec(1:int:end),PSNR_curve(6,1:int:end),'-','LineWidth',lw,'Color',colorscheme(9,:));
plot(HS_spec(1:int:end),PSNR_curve(12,1:int:end),'-','LineWidth',lw,'Color',colorscheme(6,:));
grid on
% ---- add legend
xlabel('wavelength (nm)'); ylabel('PSNR (dB)');
set(gca,'FontSize',30);
legend('naive interpolation','CNMF','FUMI','SupResPALM','HiBCD (plain CoSMF)','HiBCD (NNC CoSMF)');
axis([350 1010 20 40]);