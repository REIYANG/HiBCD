function S_Proj = Simplex_Proj(S)

S_Proj = max(bsxfun(@minus,S,max(bsxfun(@rdivide,cumsum(sort(S,1,'descend'),1)-1,(1:size(S,1))'),[],1)),0);

