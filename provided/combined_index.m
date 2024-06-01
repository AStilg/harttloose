function [out1,out2] = combined_index(in1,in2)
% combined_index.m - translates between (n,m) and combined index
%                    ci = n * (n+1) + m
%
% Usage:
% [n,m] = combined_index(ci);
% ci = combined_index(n,m);
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

% Sanity check
if nargin == 1
    if strcmp(class(in1),'logical')
        szIn1=size(in1);
        if szIn1(2)>1
            [in1,~]=find(in1);
            in1=reshape(in1,szIn1);
        else
            in1=find(in1);
        end
        
    end
    out1 = floor(sqrt(in1));
    out2 = in1 - out1.^2 - out1;
    if nargout==1
        out1=[out1;out2];
    end
elseif nargin == 2
    out1 = in1 .* (in1 + 1) + in2;
else
    error('Bad number of input/output arguments');
end

return
