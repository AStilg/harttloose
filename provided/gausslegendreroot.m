function [x,w]=gausslegendreroot(n,ab);
% gausslegendreroot - compute the roots and weight of GL quadrature.
%
% USAGE:
% [x,w] = gausslegendreroot(n);
% [x,w] = gausslegendreroot(n,ab);
%
% n  -- number of nodes
% ab -- integration limits (2 vector)
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if nargin==1
    ab=[-1,1];
end

p=[4.158083334531387e-01     2.097392773191518e-01]; %approx inverse values

[ix0]=polyval(p,n);

x0=cos(linspace(pi-1./ix0,1./ix0,n));

x=x0;
for jj=1:2
    [L]=legendrecol(n+1,x,0,1)./spherical_harmonic_normalisation([0:n+1]',0);
    dL=((-n-1)*L(end,:)+(n+1)*x.*L(end-1,:))./(1-x.^2);
    L=L(end-1,:);
    x=x-L./dL;
end

L=legendrecol(n+1,x,0,1);
wi=2.*(1-x.^2)./((n+1).^2.*L(end,:).^2/spherical_harmonic_normalisation(n+1,0)^2);

x=diff(ab)/2.*x(:)+sum(ab)/2;

w=diff(ab)/2.*wi(:);

