function [M,N,M2,N2,M3,N3,Norm,eigenvalues]=spheroidalvwf(isProlate,n,m,c,xi,eta,phi,htype);
% spheroidalvwf - compute the "unnormalised" vector spheroidal wavefunctions
% 				  of kinds (3), (4), (1) (corresponding to VSWF: 1,2,3).
%
% usage:
%
% [M,N]=spheroidalvwf(isProlate,n,m,c,xi,eta,phi,htype)
%
% [M,N]=spheroidalvwf(isProlate,n,m,c,allcoordinates,htype)
%
% [M,N,M2,M2,M3,M3]=spheroidalvwf(isProlate,n,m,c,xi,eta,phi)
%
% n, m        -- radial and azimuthal modes equal length arrays ok.
% c           -- interfocal number (not c/k)
% xi,eta, phi -- coordinates, vectors of N points.
% htype       -- 1: incoming, 2: outgoing, 3: regular (default: 3)
%
% allcoordinates -- is [N x (3,4,5)] which can also include sqrt(xi.^2-sigma)
% and sqrt(1-eta.^2) respectively.
%
% output:
%
% M,N,M2,N2,M3,N3 - if n,m are scalar then: outputs are [N x 3]. if n,m are
%					arrays then: [length(n) x 3N].
%
% PACKAGE INFO

n=n(:).';
m=m(:).';

assert(length(n)==length(m), 'n must be same length as m.')

sigma=2*(isProlate)-1;

if nargin<8
    htype=0;
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

if ~exist('xi_f','var')
    if isProlate
        xi_f=sqrt(-1+xi).*sqrt(1+xi);
    else
        xi_f=sqrt(1+xi.^2);
    end
end

eta2_f=eta_f.^2;

if isProlate
    xieta2_f=(xi-eta).*(xi+eta);
    xieta_f=sqrt(xieta2_f);
else
    xieta2_f=(xi.^2+eta.^2);
    xieta_f=sqrt(xieta2_f);
end
xi2_f=xi_f.^2;

% determine what's needed.
% find unique abs(m)...
nmax=max(n);
parities=mod(n+m,2); %1 for odd.

[uniq_m]=unique(abs(m));

R1nm=zeros(length(m),length(xi));
dR1nm=R1nm;
R2nm=R1nm;
dR2nm=R1nm;
Snm=R1nm;
dSnm=R1nm;
dphiSnm=R1nm;
eigenvalues=zeros(1,length(m)); %won't need to rotate.
Norm=zeros(1,length(m));

for ii=1:length(uniq_m)
    
    find_um=find(abs(m)==uniq_m(ii));
    parities_for_m=unique(parities(find_um));
    for jj=1:length(parities_for_m);
        find_parity=find(parities_for_m(jj)==parities(find_um));
        %vector of modes:
        
        nd=uniq_m(ii)+parities_for_m(jj)+2*[0:size(R1nm,1)-1]';
        indx=(n(find_um(find_parity))-parities_for_m(jj)-abs(m(find_um(find_parity))))/2+1;
        
        [Snmt,dSnmt,dphiSnmt,eigenvaluest,Normt]=spheroidalS1(isProlate,nmax+mod(nmax+uniq_m(ii),2)+parities_for_m(jj),uniq_m(ii),c,eta,true); %
        Snm(find_um(find_parity),:)=Snmt(indx,:);
        dSnm(find_um(find_parity),:)=dSnmt(indx,:);
        dphiSnm(find_um(find_parity),:)=dphiSnmt(indx,:);
        eigenvalues(find_um(find_parity))=eigenvaluest(indx);
        Norm(find_um(find_parity))=Normt(indx);
        
        [R1nmt,dR1nmt]=spheroidalR1(isProlate,nmax+mod(nmax+uniq_m(ii),2)+parities_for_m(jj),uniq_m(ii),c,[xi,xi_f],true);
        
        R1nm(find_um(find_parity),:)=R1nmt(indx,:);
        dR1nm(find_um(find_parity),:)=dR1nmt(indx,:);
        
        if htype<3
            %as this part is SLOW, we will only calculate for unique values
            %of xi
            [uxi,~,ic]=unique(xi);
            [R2nmt,dR2nmt]=spheroidalR2(isProlate,nmax+mod(nmax+uniq_m(ii),2)+parities_for_m(jj),uniq_m(ii),c,uxi,true);
            
            R2nm(find_um(find_parity),:)=R2nmt(indx,ic);
            dR2nm(find_um(find_parity),:)=dR2nmt(indx,ic);
        end
    end
    
end

%error check.
if length(n)==1
	Snm=Snm.';
	dSnm=dSnm.';
	dphiSnm=dphiSnm.';

	R1nm=R1nm.';
	dR1nm=dR1nm.';
	if htype<3
		R2nm=R2nm.';
		dR2nm=dR2nm.';
    end
else
	n=n';
	m=m';
    
    xi=xi.';
    eta=eta.';
    phi=phi.';
    
    xi_f=xi_f.';
    xi2_f=xi2_f.';
    
    eta_f=eta_f.';
    eta2_f=eta2_f.';
    
    xieta_f=xieta_f.';
    xieta2_f=xieta2_f.';
    
    eigenvalues=eigenvalues.';
    Norm=Norm.';
end

if ~exist('xi_f','var')
    if isProlate
        xi_f=sqrt(-1+xi).*sqrt(1+xi);
    else
        xi_f=sqrt(1+xi.^2);
    end
end


dxiR1nm=(dR1nm+(sigma-3.*xi.^2).*R1nm)./(xi.*xi2_f);
if htype<3
    dxiR2nm=(dR2nm+(sigma-3.*xi.^2).*R2nm)./(xi.*xi2_f);
end
%
if or(nargout<3,htype>0)
    
    switch htype
        case 1
            Rnm=R1nm+1i*R2nm;
            dRnm=dR1nm+1i*dR2nm;
            dxiRnm=dxiR1nm+1i*dxiR2nm;
        case 2
            Rnm=R1nm-1i*R2nm;
            dRnm=dR1nm-1i*dR2nm;
            dxiRnm=dxiR1nm-1i*dxiR2nm;
        case 3
            Rnm=R1nm;
            dRnm=dR1nm;
            dxiRnm=dxiR1nm;
    end

    M=[sigma.*(1i.*m.*eta.*Snm.*Rnm)./(xi_f.*xieta_f),   -sign(m).*(1i.*xi.*dphiSnm.*Rnm)./(xieta_f),    xi_f.*(xi.*dSnm.*Rnm-eta_f.*sigma.*eta.*Snm.*dxiRnm)./xieta2_f].*repmat((sign(m)).^m.*exp(1i*m.*phi),[1,3]);
    %we consider m correct now!
    polarTerm=m.*xi.*sign(m).*dphiSnm.*Rnm.*(1./xi2_f-1./xieta2_f)./eta_f;
    polarTerm(:,abs(eta)==1)=0;

    N=[(sigma.*((xi.^2+sigma.*eta.^2.*(1+eta.^2-3.*sigma.*xi.^2))./xieta2_f.*Snm+eta.*eta_f.*dSnm)./xieta2_f.*dxiRnm - ...
        xi.*(2*sigma.*eta./xieta2_f.*eta_f.*dSnm-(sigma.*c.^2.*eta2_f+eigenvalues).*Snm).*Rnm./xieta2_f + ...
        polarTerm).*xi_f./xieta_f, ...
        ((dSnm.*(-2*xi.^2.*xi2_f./xieta2_f.*Rnm+dRnm)./xieta2_f) - ...
        sigma.*eta.*eta_f.*Snm.*(-2.*(dRnm+(sigma-3.*xi.^2).*Rnm)./xieta2_f+(sigma.*m.^2./xi2_f-c.^2.*xi2_f+eigenvalues).*Rnm)./xieta2_f + ...
        sigma.*m.*eta.*Rnm.*sign(m).*dphiSnm./xi2_f)./xieta_f       , ...
        1i*((sign(m).*dphiSnm.*(xieta2_f.*Rnm+xi.*xi2_f.*dxiRnm)+sigma.*m.*eta.*Rnm.*dSnm)./xieta2_f./xi_f)].*repmat((sign(m)).^m.*exp(1i*m.*phi),[1,3])/c;

    M2=0;
    N2=0;
    M3=0;
    N3=0;
else
    
    Rnm=R1nm+1i*R2nm;
    dRnm=dR1nm+1i*dR2nm;
    dxiRnm=dxiR1nm+1i*dxiR2nm;
    
    M=[sigma.*(1i.*m.*eta.*Snm.*Rnm)./(xi_f.*xieta_f),   -sign(m).*(1i.*xi.*dphiSnm.*Rnm)./(xieta_f),    xi_f.*(xi.*dSnm.*Rnm-eta_f.*sigma.*eta.*Snm.*dxiRnm)./xieta2_f].*repmat((sign(m)).^m.*exp(1i*m.*phi),[1,3]);
    %we consider m correct now!
    polarTerm=m.*xi.*sign(m).*dphiSnm.*Rnm.*(1./xi2_f-1./xieta2_f)./eta_f;
    polarTerm(:,abs(eta)==1)=0;
    
    N=[(sigma.*((xi.^2+sigma.*eta.^2.*(1+eta.^2-3.*sigma.*xi.^2))./xieta2_f.*Snm+eta.*eta_f.*dSnm)./xieta2_f.*dxiRnm - ...
        xi.*(2*sigma.*eta./xieta2_f.*eta_f.*dSnm-(sigma.*c.^2.*eta2_f+eigenvalues).*Snm).*Rnm./xieta2_f + ...
        polarTerm).*xi_f./xieta_f, ...
        ((dSnm.*(-2*xi.^2.*xi2_f./xieta2_f.*Rnm+dRnm)./xieta2_f) - ...
        sigma.*eta.*eta_f.*Snm.*(-2.*(dRnm+(sigma-3.*xi.^2).*Rnm)./xieta2_f+(sigma.*m.^2./xi2_f-c.^2.*xi2_f+eigenvalues).*Rnm)./xieta2_f + ...
        sigma.*m.*eta.*Rnm.*sign(m).*dphiSnm./xi2_f)./xieta_f       , ...
        1i*((sign(m).*dphiSnm.*(xieta2_f.*Rnm+xi.*xi2_f.*dxiRnm)+sigma.*m.*eta.*Rnm.*dSnm)./xieta2_f./xi_f)].*repmat((sign(m)).^m.*exp(1i*m.*phi),[1,3])/c;
   
    Rnm=R1nm-1i*R2nm;
    dRnm=dR1nm-1i*dR2nm;
    dxiRnm=dxiR1nm-1i*dxiR2nm;
    
    M2=[sigma.*(1i.*m.*eta.*Snm.*Rnm)./(xi_f.*xieta_f),   -sign(m).*(1i.*xi.*dphiSnm.*Rnm)./(xieta_f),    xi_f.*(xi.*dSnm.*Rnm-eta_f.*sigma.*eta.*Snm.*dxiRnm)./xieta2_f].*repmat((sign(m)).^m.*exp(1i*m.*phi),[1,3]);
    %we consider m correct now!
    polarTerm=m.*xi.*sign(m).*dphiSnm.*Rnm.*(1./xi2_f-1./xieta2_f)./eta_f;
    polarTerm(:,abs(eta)==1)=0;
    
    N2=[(sigma.*((xi.^2+sigma.*eta.^2.*(1+eta.^2-3.*sigma.*xi.^2))./xieta2_f.*Snm+eta.*eta_f.*dSnm)./xieta2_f.*dxiRnm - ...
        xi.*(2*sigma.*eta./xieta2_f.*eta_f.*dSnm-(sigma.*c.^2.*eta2_f+eigenvalues).*Snm).*Rnm./xieta2_f + ...
        polarTerm).*xi_f./xieta_f, ...
        ((dSnm.*(-2*xi.^2.*xi2_f./xieta2_f.*Rnm+dRnm)./xieta2_f) - ...
        sigma.*eta.*eta_f.*Snm.*(-2.*(dRnm+(sigma-3.*xi.^2).*Rnm)./xieta2_f+(sigma.*m.^2./xi2_f-c.^2.*xi2_f+eigenvalues).*Rnm)./xieta2_f + ...
        sigma.*m.*eta.*Rnm.*sign(m).*dphiSnm./xi2_f)./xieta_f       , ...
        1i*((sign(m).*dphiSnm.*(xieta2_f.*Rnm+xi.*xi2_f.*dxiRnm)+sigma.*m.*eta.*Rnm.*dSnm)./xieta2_f./xi_f)].*repmat((-sign(m)).^m.*exp(1i*m.*phi),[1,3])/c;

    M3=[sigma.*(1i.*m.*eta.*Snm.*R1nm)./(xi_f.*xieta_f),   -sign(m).*(1i.*xi.*dphiSnm.*R1nm)./(xieta_f),    xi_f.*(xi.*dSnm.*R1nm-eta_f.*sigma.*eta.*Snm.*dxiR1nm)./xieta2_f].*repmat((sign(m)).^m.*exp(1i*m.*phi),[1,3]);
    %we consider m correct now!
    polarTerm=m.*xi.*sign(m).*dphiSnm.*R1nm.*(1./xi2_f-1./xieta2_f)./eta_f;
    polarTerm(:,abs(eta)==1)=0;
    
    N3=[(sigma.*((xi.^2+sigma.*eta.^2.*(1+eta.^2-3.*sigma.*xi.^2))./xieta2_f.*Snm+eta.*eta_f.*dSnm)./xieta2_f.*dxiR1nm - ...
        xi.*(2*sigma.*eta./xieta2_f.*eta_f.*dSnm-(sigma.*c.^2.*eta2_f+eigenvalues).*Snm).*R1nm./xieta2_f + ...
        polarTerm).*xi_f./xieta_f, ...
        ((dSnm.*(-2*xi.^2.*xi2_f./xieta2_f.*R1nm+dR1nm)./xieta2_f) - ...
        sigma.*eta.*eta_f.*Snm.*(-2.*(dR1nm+(sigma-3.*xi.^2).*R1nm)./xieta2_f+(sigma.*m.^2./xi2_f-c.^2.*xi2_f+eigenvalues).*R1nm)./xieta2_f + ...
        sigma.*m.*eta.*R1nm.*sign(m).*dphiSnm./xi2_f)./xieta_f       , ...
        1i*((sign(m).*dphiSnm.*(xieta2_f.*R1nm+xi.*xi2_f.*dxiR1nm)+sigma.*m.*eta.*R1nm.*dSnm)./xieta2_f./xi_f)].*repmat((sign(m)).^m.*exp(1i*m.*phi),[1,3])/c;
    
end
