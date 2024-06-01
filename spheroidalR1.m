function [Rnm,dRnm_mod,eigenvalues]=spheroidalR1(isProlate,n,m,c,xi_in,allvalues);
% spheroidalR1 - calculates the radial spheroidal
% 				 wavefunctions, R_{nm}(c,\eta). Can output all 
%                wavefunctions of same m and partiy.
%
% usage:
% [Rnm,dRnm_mod,eigenvalues]=spheroidalR1(isProlate,n,m,c,xi_in)
%
% or
%
% [Rnm,dRnm_mod,eigenvalues]=spheroidalR1(isProlate,n,m,c,xi_in,allvalues)
%
% Note: xi must be [N x (1 or 2)] of either: [xi(:)] or [xi(:),xifac(:)] where
% xifac is sqrt(xi.^2-(2*isProlate-1)). xifac only has useful application to prolate
% cases.
%
% Default behaviour is that only the requested value is given.
%
% NOTE: The derivative is,
%
% dRnm_mod = \frac{d}{d\xi}(\xi(\xi^2\mp 1)R_{nm}(c,\eta))
%
% NOTE: Two different series are used to generate prolate and oblate
% spheroidalR1. Though both series are valid one is more convergent than
% the other for either prolate or oblate coordinates.
%
% NOTE: if allvalues, the last few rows are guaranteed to be wrong.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if nargin<6
    allvalues=false;
end

if lower(allvalues)=='all'
    allvalues=true;
end

m=abs(m);
if n<m
    error('n must be greater or equal to |m|.')
end

% derived parameters:

xi=xi_in(:,1).';
isOdd=mod(n+abs(m),2);
N=ceil(n+abs(c)+15);

% compute superposition weights:

[U,r,eigenvalues]=spheroidal_u_coefficients(isProlate,isOdd,abs(m),c,N);

xi_arg=sqrt((1+xi).*(-1+xi));
if size(xi_in,2)==2
    xi_arg=xi_in(:,2).';
end
if ~isProlate
    xi_arg=xi;
end

if ~allvalues
    U=(U((r==n),:));
    eigenvalues=eigenvalues((r==n));
else
    n=r';
%     if (m==0)
%         U(:,1)=0;
%     end
end

%compute all functions (for derivatives too)
nv=[abs(m):r(end)+1].';
J=sbesselj(nv,c*xi_arg).';
nv=nv(r-abs(m)+1);

% compute normalisations for the bessel functions requested.
if isProlate % eta=0 for prolate
    L=legendrecol(r(end)+1,0,abs(m),1);
    besselNormalisation=L(r-abs(m)+1);
    if isOdd
        besselNormalisation=sqrt(((1-m+nv).*(1+m+nv).*(1+2.*nv))./(3+2.*nv)).*L(r-abs(m)+2);
    end
else % eta=1 for oblate:
    
    besselNormalisation=(2*r+1)./spherical_harmonic_normalisation(r,m)/sqrt(2*pi)/2;

end
spheroidalNormalisation=U*besselNormalisation;

RJs=(besselNormalisation).*J(r-abs(m)+1,:);

%compute Rnm
if isProlate
    Rnm=(xi./xi_arg).^(isOdd).*(((1i).^(r'-n').*U)*RJs)./spheroidalNormalisation;
    if m==0
        if isOdd
            Rnm(:,xi==1)=repmat((1i).^(n'+1).*(spherical_harmonic_normalisation(1,0)./spheroidalNormalisation).*U(:,1),[1,sum(xi==1)]); %hack
        else
            Rnm(:,xi==1)=repmat((1i).^(n').*(spherical_harmonic_normalisation(0,0)./spheroidalNormalisation).*U(:,1),[1,sum(xi==1)]); %hack
        end
    else
        Rnm(:,xi==1)=0;
    end
else
    Rnm=(((xi.^2)+1)./(xi.^2)).^(abs(m)/2).*(((1i).^(r'-n').*U)*RJs)./spheroidalNormalisation;

    if isOdd
        Rnm(:,xi==0)=0;
    else
        Rnm(:,xi==0)=repmat(c.^m*(((1i).^(n'-m).*U(:,1))*((besselNormalisation(1)).*(2.^(1+nv(1)).*factorial(nv(1)+1)./factorial(2*(nv(1)+1)))))./spheroidalNormalisation,[1,sum(xi==0)]);
    end
    
end

%compute dRnm_mod
if isProlate
    xifactor=(c.*xi.^2.*xi_arg);
    xifactord=(-(1+isOdd)+(3+r).*xi.^2);
    
    RJs2=besselNormalisation.*(xifactord.*J(r-abs(m)+1,:)-xifactor.*J(r-abs(m)+2,:));

    dRnm_mod=(xi./xi_arg).^(isOdd).*(((1i).^(r'-n').*U)*RJs2)./spheroidalNormalisation;

    if m==0
        if isOdd
            dRnm_mod(:,xi==1)=2*Rnm(:,xi==1);
        end
    else
        dRnm_mod(:,xi==1)=0;
    end
else
    xifactor=(xi.*(1+xi.^2));
    xifactord=((3.*xi.^2-(m-1)));
    RJs2=(besselNormalisation).*((xifactor.*r./xi+xifactord).*J(r-abs(m)+1,:)-c*xifactor.*J(r-abs(m)+2,:));
    dRnm_mod=(((xi.^2)+1)./(xi.^2)).^(abs(m)/2).*(((1i).^(r'-n').*U)*RJs2)./spheroidalNormalisation;
    
    dRnm_mod(:,xi==0)=(1+isOdd)*Rnm(:,xi==0);
end

if ~allvalues
    Rnm=Rnm(:);
    dRnm_mod=dRnm_mod(:);
end

