function imagels(x, paramstruct);
%IMAGELS is a revised version of image function in order to view the
%relative values of a data matrix.
%
%We arrange the minimum value of the matrix as 0, and the maximum value as
%100. And all the other values are linear transformation from the original
%value. Under the same color map, we can also view a local area.
%
%Inputs:
%
%    x         input matrix
%
%   paramstruct:
%         a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%    fields             Value
%
%    rowinit            initial row to draw the image
%
%    rowend             last row to draw the image
%
%    colinit            initial column to draw the image
%
%    colend             last column to draw the image
%
%(c) Copyright 2005 Lingsong Zhang (LSZHANG@email.unc.edu)

[rown, coln]=size(x);
rowinit=1;
rowend=rown;
colinit=1;
colend=coln;

if nargin>1;
    if isfield(paramstruct, 'rowinit');
        rowinit=getfield(paramstruct, 'rowinit');
    end;
    
    if isfield(paramstruct, 'rowend');
        rowend=getfield(paramstruct, 'rowend');
    end;

    if isfield(paramstruct, 'colinit');
        colinit=getfield(paramstruct, 'colinit');
    end;
    
    if isfield(paramstruct, 'colend');
        colend=getfield(paramstruct, 'colend');
    end;
end;

xmax=max(max(x));
xmin=min(min(x));

xtrans=(x-xmin).*100./(xmax-xmin);


image(xtrans(rowinit:rowend, colinit:colend));
