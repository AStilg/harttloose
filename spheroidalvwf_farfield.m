function [Mfar,Nfar,Norm]=spheroidalvwf_farfield(isProlate,n,m,c,xi,eta,phi,htype);
% spheroidalvwf_farfield - compute the far-field "unnormalised" vector spheroidal 
% 					wavefunctions of kinds (3), (4) (corresponding to VSWF: 1,2).
%
% usage:
%
% [M,N]=spheroidalvwf_farfield(isProlate,n,m,c,xi,eta,phi,htype)
%
% [M,N]=spheroidalvwf_farfield(isProlate,n,m,c,allcoordinates,htype)
%
% n, m        -- radial and azimuthal modes equal length arrays ok.
% c           -- interfocal number (not c/k)
% xi,eta, phi -- coordinates, vectors of N points.
% htype       -- 1: incoming, 2: outgoing (1 default)
%
% allcoordinates -- is [N x 3] or xi, eta, phi.
%
% if xi is infinite (or effectively so) then only the angular parts of the 
% field are reported, normalised so that e(ikr)/kr=1.
%
% is about equal accuracy to spheroidalvwf at xi=1e8. Likely this is better
% further away.
%
% output:
%
% M,N - if n,m are scalar then: outputs are [N x 3]. if n,m are
%					arrays then: [length(n) x 3N].
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

n=n(:).';
m=m(:).';

assert(length(n)==length(m), 'n must be same length as m.')

sigma=2*(isProlate)-1;

if nargin<8
    htype=1;
end

if nargin==6
    htype=eta;
end

if nargin<=6
    if size(xi,2)==3
        phi=xi(:,3);
        eta=xi(:,2);
        xi=xi(:,1);
    elseif size(xi,2)==4
        xi_f=xi(:,4);
        phi=xi(:,3);
        eta=xi(:,2);
        xi=xi(:,1);
    elseif size(xi,2)==5
        eta_f=xi(:,5);
        xi_f=xi(:,4);
        phi=xi(:,3);
        eta=xi(:,2);
        xi=xi(:,1);
    end
end

xi=xi(:);
xi(xi<isProlate+1e-15)=isProlate+1e-15; %no divide by zero hack.
eta=eta(:);
phi=phi(:);

assert(and(length(xi)==length(eta),length(eta)==length(phi)), 'Coordinates: xi, eta, phi must be same size.')

if ~exist('eta_f','var')
    eta_f=sqrt((1-eta).*(1+eta));
end


kr=c.*xi;

if kr<1e5
    warning('Small values of kr are known to be inaccurate, relative error is proportional to 1/sqrt(kr)!');
end

sgni=2*htype-3;
eikronkr=exp(-sgni*1i*kr)./kr;

eikronkr(kr>1e300)=1;

eta2_f=eta_f.^2;


nmax=max(n);
parities=mod(n+m,2); %1 for odd.

[uniq_m]=unique(abs(m));

Snm=zeros(length(m),length(xi));
dSnm=Snm;
dphiSnm=Snm;
eigenvalues=zeros(1,length(m)); %won't need to rotate.
Norm=zeros(1,length(m));

for ii=1:length(uniq_m)
    
    find_um=find(abs(m)==uniq_m(ii));
    parities_for_m=unique(parities(find_um));
    for jj=1:length(parities_for_m);
        find_parity=find(parities_for_m(jj)==parities(find_um));
        %vector of modes:
        
        nd=uniq_m(ii)+parities_for_m(jj)+2*[0:size(Snm,1)-1]';
        indx=(n(find_um(find_parity))-parities_for_m(jj)-abs(m(find_um(find_parity))))/2+1;
        
        [Snmt,dSnmt,dphiSnmt,eigenvaluest,Normt]=spheroidalS1(isProlate,nmax+mod(nmax+uniq_m(ii),2)+parities_for_m(jj),uniq_m(ii),c,eta,true); %
        Snm(find_um(find_parity),:)=Snmt(indx,:);
        dSnm(find_um(find_parity),:)=dSnmt(indx,:);
        dphiSnm(find_um(find_parity),:)=dphiSnmt(indx,:);
        eigenvalues(find_um(find_parity))=eigenvaluest(indx);
        Norm(find_um(find_parity))=Normt(indx);
        
    end
    
end


% Mfar=(sgni).^(n+1)*[0*dSnm,1i*sign(m)*dphiSnm,-dSnm].*(-sign(m)).^m.*exp(1i*m*phi)
% Nfar=(sgni).^(n  )*[0*dSnm,dSnm,1i*sign(m)*dphiSnm].*(-sign(m)).^m.*exp(1i*m*phi)

m=m';
n=n';
eigenvalues=eigenvalues.';
xi=xi.';
eta=eta.';
phi=phi.';
eta2_f=eta2_f.';

Mfar=((1i*sgni).^(n+1).*repmat(eikronkr.'.*(sign(m)).^m.*exp(1i*m.*phi),[1,3])).*[sigma.*(1i.*m.*eta.*Snm)./xi.^2,-1i*sign(m).*dphiSnm,dSnm];
Nfar=((1i*sgni).^(n  ).*repmat(eikronkr.'.*(sign(m)).^m.*exp(1i*m.*phi),[1,3])).*[1i*(sigma.*c.^2.*eta2_f+eigenvalues).*Snm./xi/c,dSnm,1i*sign(m).*dphiSnm];
