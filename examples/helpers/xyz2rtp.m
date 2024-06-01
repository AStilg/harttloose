function [r,theta,phi] = xyz2rtp(x,y,z)
% xyz2rtp.m
% Coordinate transformation from cartesian to spherical
% theta is the polar angle, measured from the +z axis,
% and varies from 0 to pi
% phi is the azimuthal angle, measured from the +x axis, increasing
% towards the +y axis, varying from 0 to 2*pi
%
% Usage:
% [r,theta,phi] = xyz2rtp(x,y,z);
% where x, y, z, r, theta, phi are all scalars or
% equal length vectors
% or
% r = xyz2rtp(x);
% where x = [ x y z ] and r = [ r theta phi ]
%
% Angles are in radians
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if nargin == 1
   y = x(:,2);
   z = x(:,3);
   x = x(:,1);
end

xy = sqrt( x.*x + y.*y );
theta = mod(atan2(xy,z)+2*pi,2*pi);
phi = mod(atan2(y,x)+2*pi,2*pi);
r = sqrt( x.*x + y.*y + z.*z );

if nargout < 2
   r = r(:);
   theta = theta(:);
   phi = phi(:);
   r = [ r theta phi ];
end

return
