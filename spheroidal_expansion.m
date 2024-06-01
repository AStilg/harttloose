function [U,Nnm0,Nnmc,im_v]=spheroidal_expansion(isProlate,c,nmax);
% spheroidal_expansion - outputs all mode couplings between spherical and
% 							spheroidal funcions up to nmax for
%							the electromagnetic problem.
%
% usage:
% [U,Nnm0,Nnmc,im_vector]=spheroidal_expansion(isProlate,c,nmax);
%
% U        -- scalar function expansion up to 2*(abs(c)+nmax)
% Nnm0     -- vector spherical normalisation
% Nnmc     -- vector spheroidal "normalisation"
% im_v 	   -- phase factor
% iProlate -- true/false
% c        -- half interfocal distance
% nmax     -- see note on output U
%
% vector operator for number of orders_wanted:
% U_v = U(1:orders_wanted,2:length(Nnm0)).*im_v(1:orders_wanted).*(im_v(2:length(Nnm0))./Nnm0(2:length(Nnm0)))';
% inverse vector operator:
% iU_v= U(1:orders_wanted,2:length(Nnm0))'.*im_v(1:orders_wanted)'.*(im_v(2:length(Nnm0)).*Nnm0(2:length(Nnm0)));
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

N=round(abs(c)+nmax);

U=sparse(nmax*(nmax+2)+1,nmax*(nmax+2)+1);

for isOdd=0:1
    [u,n]=spheroidal_u_coefficients(isProlate,isOdd,0,c,ceil((N)/2));
	
    validValues=(n<=nmax);
    ci=combined_index(n(validValues),0)+1;
	
    U(ci,ci)=u(1:sum(validValues),1:sum(validValues));
	
    for m=1:N
	
        [u,n]=spheroidal_u_coefficients(isProlate,isOdd,m,c,ceil((N-m)/2));
        
		validValues=(n<=nmax);
		ci=combined_index(n(validValues),0)+1;
	
        ci=combined_index(n(n<=nmax),m)+1;
        U(ci,ci)=u(1:sum(validValues),1:sum(validValues));
        
        ci=combined_index(n(n<=nmax),-m)+1;
        U(ci,ci)=u(1:sum(validValues),1:sum(validValues));
        
    end
end

[n,~]=combined_index([0:size(U,1)-1]');

Nnm0=sqrt(n.*(n+1));
Nnmc=sqrt(sum(U.^2.*sqrt(n.*(n+1))',2));
im_v=1i.^(n);