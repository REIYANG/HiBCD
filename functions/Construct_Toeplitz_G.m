function [G,B,S] = Construct_Toeplitz_G(M,N,ks,sigVal,dsRatio)
kernel = fspecial('Gaussian',[1 ks],sigVal);
h1 = [kernel((ks+1)/2:end) zeros(1,M-(ks+1)/2)];
h2 = [kernel((ks+1)/2:end) zeros(1,N-(ks+1)/2)];
B1 = sparse(toeplitz(h1,h1)); 
B2 = sparse(toeplitz(h2,h2));
D1 = sparse(kron(eye(M/dsRatio),[1;zeros(dsRatio-1,1)]));
D2 = sparse(kron(eye(N/dsRatio),[1;zeros(dsRatio-1,1)]));
B = KernelToMatrix(kernel,M,N);
G = sparse(kron(B2*D2,B1*D1));
S = sparse(kron(D2,D1));
end

function [B]=KernelToMatrix(KerBlu,nr,nc)
% flip the kernel
KerBlu=rot90(KerBlu,2);
mid_col=round((nc+1)/2);
mid_row=round((nr+1)/2);
% the size of the kernel
[len_hor,len_ver]=size(KerBlu);
lx = (len_hor-1)/2;
ly = (len_ver-1)/2;
B=zeros(nr,nc);
% range of the pixels
B(mid_row-lx:mid_row+lx,mid_col-ly:mid_col+ly)=KerBlu;
B=circshift(B,[-mid_row+1,-mid_col+1]);
end