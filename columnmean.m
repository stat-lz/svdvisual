function y=columnmean(x);
%COLUMNMEAN column mean of data matrix
%
%  y=columnmean(x)
%
%Copyright(c) 2004 Lingsong Zhang(LSZHANG@email.unc.edu)
[nrow, ncol]=size(x);
columnmeanvec=mean(x);
onetemp=ones(nrow, 1);
y=onetemp*columnmeanvec;