function [Rnm,dRnm_mod,eigenvalues,Rnm_full,sumz]=spheroidalR2(isProlate,n,m,c,xi,allvalues);
% spheroidalR2 - calculates the radial spheroidal
% wavefunctions, R_{nm}(c,\eta). Can output all wavefunctions of
% same m and partiy.
%
% usage:
% [Rnm,dRnm_mod,eigenvalues]=spheroidalR2(isProlate,n,m,c,xi)
%
% or
%
% [Rnm,dRnm_mod,eigenvalues]=spheroidalR2(isProlate,n,m,c,xi,allvalues)
%
% valid input arguments:
% prolate: xi>=1.
% oblate: xi>=0.
%
% Default behaviour is that only the requested value is given.
%
% NOTE: The derivative is,
%
% dRnm_mod = \frac{d}{d\xi}(\xi(\xi^2-sign(c.^2))R_{nm}(c,\eta))
% thus:
% dR=(dRnm_mod+(sign(c.^2)-3.*xi.^2).*Rnm)./(xi.*(xi.^2-sign(c.^2)))
%
% NOTE: if allvalues, the last few rows are guaranteed to be wrong.
%
% PACKAGE INFO

if nargin<6
    allvalues=false;
end

if lower(allvalues)=='all'
    allvalues=true;
end

if n<abs(m)
    error('n must be greater or equal to |m|.')
end
m=abs(m);

% derived parameters:
xi=(xi(:).');
isOdd=mod(n+abs(m),2); %not needed for this method
N=ceil(n+abs(c)+15);

% compute superposition weights:
[U,r,eigenvalues]=spheroidal_u_coefficients(isProlate,isOdd,abs(m),c,N);

% Ifac=U0(r==n,r==nd);

%use clenshaw-curtis
[acoseta,wi]=gausslegendreroot(4*(N+1),[0,pi/2]); %it'll be fine.
eta=cos(acoseta);

[XI,ETA]=meshgrid(xi,eta);

%convert coordinates to oblate is needed.
sigma=2*isProlate-1;

%compute all functions (for derivatives too)

nd=r';
% nd(nd>n)=[]; %cut

n_for_j=unique([nd,nd+1]);
indx=nd-n_for_j(1)+1;

Jnd=sbessely(n_for_j,kr(XI,ETA,c,sigma)); %all eta xi n
dJnd=c.^2.*XI(:).*(nd.*Jnd(:,indx)./(kr(XI(:),ETA(:),c,sigma)).^2-Jnd(:,indx+1)./(kr(XI(:),ETA(:),c,sigma))); %all derivatives
Pnd=legendrecol(max(nd(:))+1,costh(XI,ETA,c,sigma),m,1);

dPnd=sigma.*c.^2.*ETA(:).*((-1-nd).*ETA(:).*XI(:).*Pnd(nd-abs(m)+1,:).'+sqrt((1-m+nd).*(nd+m+1).*(2*nd+1)./(2*nd+3)).*kr(XI(:),ETA(:),c,sigma).*Pnd(nd-abs(m)+2,:).'./c)./(kr(XI(:),ETA(:),c,sigma)).^2./(-sigma+XI(:).^2);

PJ1=(Jnd(:,indx).*Pnd(nd-abs(m)+1,:)');
dPJ1=(3*(XI(:)).^2-sigma).*PJ1 + (XI(:)).*(XI(:).^2-sigma).*(Pnd(nd-abs(m)+1,:).'.*dJnd+Jnd(:,indx).*dPnd);
% dPJ1=reshape(Pnd(end-1,:).'.*dJnd+Jnd(:,1).*dPnd,size(XI));

temp=legendrecol(r(end)+1+2,eta,abs(m),1);
Y=2*temp(r-abs(m)+1,:).*(sin(acoseta).*wi)'; %may as well apply weights too

PJ1=reshape(PJ1,[size(XI),size(PJ1,2)]);
dPJ1=reshape(dPJ1,[size(XI),size(dPJ1,2)]);

% integrals=zeros([size(U,1),size(PJ1,2),length(xi)]);
% dintegrals=integrals;

MATINTS=zeros(size(Y,1),size(PJ1,2),size(PJ1,3));
dMATINTS=MATINTS;
for ii=1:size(PJ1,3)
    MATINTS(:,:,ii)=Y*PJ1(:,:,ii);
    dMATINTS(:,:,ii)=Y*dPJ1(:,:,ii);
end

sums=zeros(size(U,1),size(MATINTS,2),size(MATINTS,3));
dsums=sums;
for ii=1:size(MATINTS,3)
    U2=U;
    
%     if (m==0)
%         U2(:,1)=0;
%     end
    
    sums(:,:,ii)=U2*MATINTS(:,:,ii);
    dsums(:,:,ii)=U2*dMATINTS(:,:,ii);
end

sumz=permute(sums,[1,3,2]);
dsumz=permute(dsums,[1,3,2]);


I2=(1i).^(r'-r).*U/2/pi;

Rnm=zeros(length(r),length(xi));
dRnm_mod=Rnm;
nref=[1:length(r)]';

for ii=1:length(xi)
    
    Rnm_full=sumz(:,:,ii)./I2;
    dRnm_mod_full=dsumz(:,:,ii)./I2; %need to selah
    
    est_offset=floor((r-m)./(8.5*(abs(xi(ii))-isProlate)+1)/2);
    col=nref-est_offset;
    col(col<1)=1;
    
    R_t=tril(Rnm_full(:,:));
    dR_t=tril(dRnm_mod_full(:,:));

    Rnm(:,ii)=R_t(sub2ind(size(R_t),nref,col));
    dRnm_mod(:,ii)=dR_t(sub2ind(size(dR_t),nref,col));
end
if ~allvalues
    Rnm=Rnm(r==n,:);
    dRnm_mod=dRnm_mod(r==n,:);
end
end

function kr=kr(xi,eta,c,sigma);
kr=c.*sqrt(xi.^2+sigma.*(eta.^2-1));
end

function costh=costh(xi,eta,c,sigma);
costh=eta.*xi./sqrt(xi.^2+sigma.*(eta.^2-1));
end