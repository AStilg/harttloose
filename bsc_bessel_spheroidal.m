function [nn,mm,a,b] = bsc_bessel_spheroidal(isProlate,c,nmax,theta,lmode,Etheta,Ephi);
% bsc_bessel_spheroidal - Finds VSWF representation of bessel wave, output BSCs producing a 
% 						  beam with amplitude determined by the field input. 
% 						  If a vector input then a and b are matrices.
%
% Usage:
% [n,m,a,b] = bsc_bessel_spheroidal(isProlate,c,nmax,theta,lmode,Etheta,Ephi);
%
% or
%
% [a,b] = bsc_bessel_spheroidal(isProlate,c,nmax,theta,lmode,Etheta,Ephi);
%
% isProlate - true/false
% c         - interfocal parameter that includes (lambda normalised) wavenumber
% theta - gives the cone angle (in radians) of the wave.
% lmode - azimuthal mode
% Etheta - polarisation in +x
% Ephi   - polarisation in +y
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

%special case ... linear polarisation

if or(size(theta,2)~=1,size(lmode,1)~=size(theta,1))
    error('ott:bsc_bessel:badinput','Input is not valid. Inputs must be column vectors.')
end

if nargin==5
    %this is sufficient.
    Etheta=lmode;
    lmode=0*theta;
end

if nargin==6
    %two cases distinguishing only is ETheta is NOT wide.
    Ephi=0;
    if size(Etheta,2)==1
        Ephi=Etheta;
        Etheta=lmode;
        lmode=0*theta;
    else
        %do nothing
    end
end

%nargin==7 is default behaviour

%% preamble
nTheta=length(theta);

Etheta=Etheta.*(1i).^(lmode+1);%.*sign(cos(theta(indexes)))
Ephi=Ephi.*(1i).^(lmode+1);


if any(theta>pi)
    warning('ott:bsc_plane:thetatoobig','Theta is larger than PI, treating all elements as angles in degrees.')
    theta=theta/180*pi;
end

tetms=0;
if size(Etheta,2)==2
    % lmode=lmode+[-1,1]; %left is -1i for Ephi, right is +1i for Ephi;
    % lmode=lmode(:);
    tetms=[-1,1];
    indexes=[[1:size(theta,1)].';[1:size(theta,1)].']; %packing for theta uniform
    
    %need to record polarisation and modify Etheta and Ephi.
    polarisation_weights=([1,-1i;1,1i]*(Etheta.')).'/2;
    Ephi=[-1i*ones(size(theta,1),1);1i*ones(size(theta,1),1)].*polarisation_weights(:);
    Etheta=[ones(size(theta,1),1);ones(size(theta,1),1)].*polarisation_weights(:).*sign(cos(theta(indexes)));%.*sign(cos(theta(indexes)))
    
end

nvec=0:nmax(1);
if size(nmax,2)==2
    nmax=nmax(1);
    nvec=nmax;
end

orders_wanted=nmax*(nmax+2)+1;

nmax_int=floor(nmax+abs(c)+10);
[u,Nnm0,~,im_v]=spheroidal_expansion(isProlate,c,nmax_int);

orders_summed=nmax_int*(nmax_int+2)+1;
U = u(1:orders_wanted,2:orders_summed)./Nnm0(2:orders_summed)'.^2; 

%calculate the mode indices we are going to find.
[nn,mm]=combined_index([0:orders_wanted-1]');

DTY=zeros(2*nmax_int*(nmax_int+2),size(theta,1));
DPY=DTY;

a=sparse([size(U,1)],size(theta,1));
b=sparse([size(U,1)],size(theta,1));

eta=cos(theta);
sizetheta=length(theta);
ulmode=unique(lmode);
for jj=1:length(ulmode);
    m0=ulmode(jj);
    populated_cols=find(lmode==m0);
    for ii=1:length(tetms)
        m=m0+tetms(ii);
        DTY=0*DTY;
        DPY=0*DPY;
        for sigma=0:1
            [~,dtS,dpS]=spheroidalS1(isProlate,nmax_int+mod(nmax_int+m,2)+sigma,abs(m),0,eta,true);
            n=2*[0:size(dtS,1)-1]'+sigma+abs(m);
            keep=(n<=floor(sqrt(size(a,1))));

            DTY(combined_index(n(keep),m)+1,:)=-dtS((keep),:).*(sign(m)).^(m);
            DPY(combined_index(n(keep),m)+1,:)=sign((m))*1i*dpS((keep),:).*(sign(m)).^(m);
        end
        DTS=(U*(DTY([1:size(U,2)]+1,:)));
        DPS=(U*(DPY([1:size(U,2)]+1,:)));

        a(:,populated_cols)=a(:,populated_cols)+4*pi*(-1i).^(nn+1).*(conj(DPS(:,populated_cols)).*Etheta((ii-1)*sizetheta+populated_cols).'-conj(DTS(:,populated_cols)).*Ephi((ii-1)*sizetheta+populated_cols).');
        b(:,populated_cols)=b(:,populated_cols)+4*pi*(-1i).^(nn)  .*(+conj(DTS(:,populated_cols)).*Etheta((ii-1)*sizetheta+populated_cols).'+conj(DPS(:,populated_cols)).*Ephi((ii-1)*sizetheta+populated_cols).');
    end
end
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
