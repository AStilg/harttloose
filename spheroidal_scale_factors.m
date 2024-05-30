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
% PACKAGE INFO

sigma=2*isProlate-1;

h=conk*[sqrt(xi(:).^2-eta(:).^2.*sigma)./sqrt(xi(:).^2-sigma),sqrt(xi(:).^2-eta(:).^2.*sigma)./sqrt(1-eta(:).^2),sqrt((1-eta(:).^2).*(xi.^2-sigma))];