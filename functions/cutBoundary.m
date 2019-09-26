function Y_cut = cutBoundary(Y,L1,L2,kernelSize)

M = size(Y,1);
img = reshape(Y',L1,L2,[]);
margin = (kernelSize-1)/2;
img_cut = img(margin+1:end-margin,margin+1:end-margin,:);
Y_cut = reshape(img_cut,[],M)';