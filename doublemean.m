function y=doublemean(x);
%DOUBLEMEAN Lingsong's Double Mean for SVD decomposition
%
%  doublemean(x)=rowmean(x)+colmean(x)-overallmean(x);
%
% return a matrix with entry listed above.
%
%Copyright(c) 2004, Lingsong Zhang(LSZHANG@email.unc.edu)

y=rowmean(x)+columnmean(x)-overallmean(x);
