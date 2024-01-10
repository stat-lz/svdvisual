function y=rowmean(x);
%ROWMEAN row mean of data matrix
%
%  y=rowmean(x);
%
%Copyright(c) 2004 Lingsong Zhang(LSZHANG@email.unc.edu)
[nrow, ncol]=size(x);
rowmeanvec=mean(x');
onetemp=ones(1, ncol);
y=rowmeanvec'*onetemp;