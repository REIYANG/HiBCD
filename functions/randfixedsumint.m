function R = randfixedsumint(m,n,S)
% This generates an m by n array R.  Each row will sum to S, and
% all elements are all non-negative integers.  The probabilities
% of each possible set of row elements are all equal.
% RAS - Mar. 4, 2017
if ceil(m)~=m|ceil(n)~=n|ceil(S)~=S|m<1|n<1|S<0
    error('Improper arguments')
else
    P = ones(S+1,n);
    for in = n-1:-1:1
        P(:,in) = cumsum(P(:,in+1));
    end
    R = zeros(m,n);
    for im = 1:m
        s = S;
        for in = 1:n
            R(im,in) = sum(P(s+1,in)*rand<=P(1:s,in));
            s = s-R(im,in);
        end
    end
end