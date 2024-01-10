function [u, s, v]=svdls(data, nosvd)
%SVD decomposition, a revised version of the function svds.m, which is
%defined in the MATLAB system. 
%
%The revise version is used to set the u, s, v matrix to be unique by some
%constraint. See Zhang (2006) for details.
%
%  [u, s, v]=svdls(data, nosvd);
%
%(c)Copyright Lingsong Zhang (lszhang@unc.edu) 2006

datarank=rank(data);
if nargin==1;
   nosvd=datarank;
elseif nosvd>datarank;
    nosvd=datarank;
end; %check whether the nosvd is specified or mis-specified.

[u, s, v]=svds(data, nosvd);

for isvd=1:nosvd;
    vtemp=v(:, isvd);
    [vsort, vindex]=max(abs(vtemp));
    vindex=min(vindex);
    u(:, isvd)=sign(vtemp(vindex)).*u(:, isvd);
    v(:, isvd)=sign(vtemp(vindex)).*v(:, isvd);
end;




