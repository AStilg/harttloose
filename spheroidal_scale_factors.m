function h = spheroidal_scale_factors(isProlate,conk,xi,eta,phi);
% spheroidal_scale_factors - creates the scale factors needed for vector
% 							operations of spheroidal vector functions and 
%							integrals.
%
% usage:
%
% h = spheroidal_scale_factors(isProlate,conk,xi,eta,phi)
%
% h - spheroidal coordinate weights for integration, vector operations.
% isProlate - prolate (true) oblate (false)
% conk - interfocal distance.
% xi - pseudo-radial coordinate
% eta - pseudo-angular coordinate
% phi - azimuthal angle (is not used in calculation)
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

sigma=2*isProlate-1;

h=conk*[sqrt(xi(:).^2-eta(:).^2.*sigma)./sqrt(xi(:).^2-sigma),sqrt(xi(:).^2-eta(:).^2.*sigma)./sqrt(1-eta(:).^2),sqrt((1-eta(:).^2).*(xi.^2-sigma))];