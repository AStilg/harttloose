function Nmax = ka2nmax(ka)
% ka2nmax.m - Finds a reasonable maximum order to truncate at given
%             a size parameter ka
%
% Returns Nmax = ka + 3 (ka)^(1/3)
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

Nmax = ka + 3 * ka.^(1/3);
Nmax = ceil(Nmax);

return
