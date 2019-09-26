load('Chikusei.mat'); % the data file downloaded from http://naotoyokoya.com/Download.html

Lx = 1080;
Ly = 1080;
chikusei = chikusei(480:Lx+480,580:Ly+580,:);
chikusei = chikusei-min(chikusei(:));
chikusei = chikusei/max(chikusei(:));

save('Chikusei_1080.mat','chikusei');