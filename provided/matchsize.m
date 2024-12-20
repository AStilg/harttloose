function [A,B,C] = matchsize(A,B,C)
% matchsize.m - checks that all vector inputs have the same
%     number of rows, and expands single-row inputs by repetition
%     to match the input row number.
%
% Usage:
% [A,B] = matchsize(A,B)
% [A,B,C] = matchsize(A,B,C)
%
% Either 2 or 3 input vectors/scalars are allowed.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

An = size(A);
Bn = size(B);
if nargin >= 3
   Cn = size(C);
else
   Cn = [1 1];
   C = [0];
end

nmax = max(An(1),Bn(1));
nmax = max(nmax,Cn(1));

if An(1) < nmax
   if An(1) == 1
      A = ones(nmax,1) * A;
  else
      error(['Number of rows in inputs must be one or equal.']);
  end
end
if Bn(1) < nmax
   if Bn(1) == 1
      B = ones(nmax,1) * B;
  else
      error(['Number of rows in inputs must be one or equal.']);
  end
end
if Cn(1) < nmax & nargin >= 3
   if Cn(1) == 1
      C = ones(nmax,1) * C;
  else
      error(['Number of rows in inputs must be one or equal.']);
  end
end

return
