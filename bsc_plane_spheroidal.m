function [nn,mm,a,b] = bsc_plane_spheroidal(isProlate,c,nmax,theta,phi,Etheta,Ephi)
% bsc_plane_spheroidal - Finds VSWF representation of plane wave, output 
% 						 BSCs producing a plane with amplitude determined 
% 						 by the field input. If a vector input then a and 
%						 b are matrices.
%
% Usage:
% [n,m,a,b] = bsc_plane_spheroidal(isProlate,c,nmax,theta,phi,Etheta,Ephi)
%
% or
%
% [a,b] = bsc_plane_spheroidal(isProlate,c,nmax,theta,phi,Etheta,Ephi)
%
% theta and phi (in radians) give the direction of propagation
% of the plane wave. +z direction is theta = 0, phi = any value
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

%special case ... linear polarisation
if nargin==6
    szE=size(Etheta);
    if szE(2)==2
        if szE(1)==1
            Etheta=repmat(Etheta,[length(phi),1]);
        end
        Ephi=-sin(phi).*Etheta(:,1)+cos(phi).*Etheta(:,2);
        Etheta=sign(cos(real(theta))).*(cos(phi).*Etheta(:,1)+sin(phi).*Etheta(:,2)); %sign(cos(real(theta))).*
    else
        error('hartm:bsc_plane_spheroidal:badargs','Most likely: Invalid uniform polarisation. Otherwise: Invalid argument.');
    end
end

theta=theta(:);
phi=phi(:);
Etheta=1i*(Etheta(:).');
Ephi=1i*(Ephi(:).');

if any(theta>pi)
    warning('hartm:bsc_plane_spheroidal:thetatoobig','Theta is larger than PI, treating all elements as angles in degrees.')
    theta=theta/180*pi;
end

orders_wanted=nmax*(nmax+2)+1;

nmax_int=floor(nmax+abs(c)+10);
orders_summed=nmax_int*(nmax_int+2)+1;
[u,Nnm0,~,im_v]=spheroidal_expansion(isProlate,c,nmax_int);

U = u(1:orders_wanted,2:orders_summed)./Nnm0(2:orders_summed)'.^2; 

%calculate the mode indices we are going to find.
[nn,mm]=combined_index([0:orders_wanted-1]');

DTY=zeros(2*nmax_int*(nmax_int+2),length(theta));
DPY=DTY;

eta=cos(theta);
for sigma=0:1
    m=0;
    
    [~,dtS,~]=spheroidalS1(isProlate,nmax_int+mod(nmax_int+m,2)+sigma,abs(m),0,eta,true);
    n=2*[0:size(dtS,1)-1]'+sigma+abs(m);
    DTY(combined_index(n,m)+1,:)=-dtS;
    
    for m = 1:nmax_int

        [~,dtS,dpS]=spheroidalS1(isProlate,nmax_int+mod(nmax_int+m,2)+sigma,abs(m),0,eta,true);
        n=2*[0:size(dtS,1)-1]'+sigma+abs(m);

        DTY(combined_index(n,m)+1,:)=-dtS;
        DPY(combined_index(n,m)+1,:)=1i*dpS;
        DTY(combined_index(n,-m)+1,:)=-dtS.*(-1).^m;
        DPY(combined_index(n,-m)+1,:)=-1i*dpS.*(-1).^m;
    
    end
end

DTS=(U*(DTY([1:size(U,2)]+1,:))).*exp(1i*mm*phi.');
DPS=(U*(DPY([1:size(U,2)]+1,:))).*exp(1i*mm*phi.');

a=4*pi*(-1i).^(nn+1).*(conj(DPS).*(Etheta)-conj(DTS).*(Ephi));
b=4*pi*(-1i).^(nn)  .*(+conj(DTS).*(Etheta)+conj(DPS).*(Ephi));

p=abs(a).^2+abs(b).^2;
binaryvector=(p>1e-15*max(p));

if nargout>2
    nn=nn(any(binaryvector,2),:);
    mm=mm(any(binaryvector,2),:);
    a=a(any(binaryvector,2),:);
    b=b(any(binaryvector,2),:);
end

if nargout==2
    a(~binaryvector)=0;
    b(~binaryvector)=0;
    nn=sparse(a);
    mm=sparse(b);
end
return
