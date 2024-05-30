function [T,T2,c,CM1,CM2,CM3]=stmatrix_spheroid_pm(nmax,k_medium,k_particle,ac);
% stmatrix_spheroid_pm - computes point-matched sT-matrix for spheroids
% using point-matching method.
%
% usage:
% [T,T2,c_d]=stmatrix_spheroid_pm(nmax,k_medium,k_particle,ac);
%
% where,
% nmax       - the number of "radial" basis functions to terminate. Found using
%              Wiscombe or Brock condition.
% k_medium   - wave number in the medium.
% k_particle - wave number in the particle.
% ac         - [semi-major,semi-minor] distance (oblate), swap for prolate.
%
% outputs,
% T    - Tmatrix.
% T2   - RgQ matrix.
% c_d  - interfocal distance (in units of ac).
%
% PACKAGE INFO

%% testing
% aspect_ratio=3;
% major_axis=1+0/64; %wavelength in medium
% minor_axis=major_axis/aspect_ratio;
% ac=[minor_axis,major_axis];
% k_medium=2*pi*1;
% k_particle=2*pi*1.2;
% nmax=round(max(ac)*k_particle+3);
%% setup
nmaxs=nmax;
if numel(nmaxs)==1
    nmaxs=[nmax,nmax];
end
nmax=max(nmaxs);

major_axis=ac(2);
minor_axis=ac(1);
aspect_ratio=major_axis/minor_axis;

k_rel=k_particle/k_medium;

isProlate=true;
if ac(2)<ac(1)
	isProlate=false;
end
sigma=2*isProlate-1;

c=(sqrt(sigma*(aspect_ratio^2-1))*minor_axis); % conk_m
cp=(k_rel*sqrt(sigma*(aspect_ratio^2-1))*minor_axis);% conk_p

% create match points
match_points=2*(2*nmax+2)+100;

% eta=linspace(-1+1/match_points,1-1/match_points,match_points);
eta=gausslegendreroot(match_points);
% eta=cumsum(rand(1,match_points)*4/match_points)-1;
eta=eta(:)';

eta(eta<0)=[];
match_points=length(eta);

xi=abs(sqrt(sigma*ac(1).^2./c.^2+1));

% create factors.
eta2_f=((1-eta).*(1+eta));
eta_f=sqrt(eta2_f);

if isProlate
    xi_f=sqrt(-1+xi).*sqrt(1+xi);
    xieta2_f=(xi-eta).*(xi+eta);
    xieta_f=sqrt(xieta2_f);
else
    xi_f=sqrt(1+xi.^2);
    xieta2_f=(xi.^2+eta.^2);
    xieta_f=sqrt(xieta2_f);
end
xi2_f=xi_f.^2;

%% pack matrix of fields:

total_orders=nmax*(nmax+2);
ci_total=[0:total_orders]';
lci=length(ci_total);

kc=k_medium*c;
kcp=k_medium*cp;
%create modes based on loops in m

[n,m]=combined_index(ci_total);

XI=xi*ones(size(eta(:)));

[M2t,N2t,~,~,M1t,N1t]=spheroidalvwf(isProlate,n,m,kc,XI,eta(:),0*eta(:));
[M3t,N3t]=spheroidalvwf(isProlate,n,m,kcp,XI,eta(:),0*eta(:),3);
M2t=M2t.';
N2t=N2t.';
M1t=M1t.';
N1t=N1t.';
M3t=M3t.';
N3t=N3t.';
for ii=1:size(M1t,2)
    M1t(:,ii)=reshape(cross([1,0,0].*ones(size(M1t,1)/3,1),reshape(M1t(:,ii),[],3)),[],1);
    N1t(:,ii)=reshape(cross([1,0,0].*ones(size(M1t,1)/3,1),reshape(N1t(:,ii),[],3)),[],1);
    M2t(:,ii)=reshape(cross([1,0,0].*ones(size(M1t,1)/3,1),reshape(M2t(:,ii),[],3)),[],1);
    N2t(:,ii)=reshape(cross([1,0,0].*ones(size(M1t,1)/3,1),reshape(N2t(:,ii),[],3)),[],1);
end

for ii=1:size(M3t,2)
    M3t(:,ii)=reshape(cross([1,0,0].*ones(size(M3t,1)/3,1),reshape(M3t(:,ii),[],3)),[],1);
    N3t(:,ii)=reshape(cross([1,0,0].*ones(size(M3t,1)/3,1),reshape(N3t(:,ii),[],3)),[],1);
end

M1t=M1t(end/3+1:end,:);
N1t=N1t(end/3+1:end,:);
M2t=M2t(end/3+1:end,:);
N2t=N2t(end/3+1:end,:);
M3t=M3t(end/3+1:end,:);
N3t=N3t(end/3+1:end,:);

CM1=[[M1t;N1t],[N1t;M1t]];
CM2=-[[M2t;N2t],[N2t;M2t]];
CM3=[[M3t;kcp/kc*N3t],[N3t;kcp/kc*M3t]];


symmetry=[inf,1];
modeMatrix=mode_couplings_spheroidal(nmaxs(1),symmetry);
modeMatrixI=mode_couplings_spheroidal(nmaxs(2),symmetry);


% imagesc(modeMatrix)

T=sparse(2*lci,2*lci,(lci)^2);
T2=T;
% warning('off');
% index=find(modeMatrix(:,1))
for ii=1:lci
    index=find(modeMatrix(:,ii));
    
    index2=find(modeMatrixI(:,lci+ii));
%     [~,m]=combined_index(ii-1);
    
%     if m==0    
%        index(index>nmax*(nmax+2)+1)=[];
%        index2(index2<nmax*(nmax+2)+1)=[];
%     end

    %if m=0
    
    Ta=[CM2(:,index),CM3(:,index)]\(CM1(:,ii)); %E fields
    Tb=[CM2(:,index2),CM3(:,index2)]\CM1(:,ii+lci); %H fields
    
    % Ta=lsqminnorm([CM2(:,index),CM3(:,index)],(CM1(:,ii)));
    % Tb=lsqminnorm([CM2(:,index2),CM3(:,index2)],(CM1(:,ii+lci)));
    
%     Ta=lsqminnorm(coeff_matrix(:,[index;index+2*lci]),(incident_wave_matrix(:,ii))); %E fields
%     Tb=lsqminnorm(coeff_matrix(:,[index2;index2+2*lci]),(incident_wave_matrix(:,ii+lci))); %H fields
%     Ta=linsolve(coeff_matrix(:,[index;index+2*lci]),(incident_wave_matrix(:,ii)));
%     Tb=linsolve(coeff_matrix(:,[index2;index2+2*lci]),(incident_wave_matrix(:,ii+lci))); %H fields

    T(index,ii)=Ta(1:end/2);
    T(index2,ii+lci)=Tb(1:end/2);
    
    T2(index,ii)=Ta(end/2+1:end);
    T2(index2,ii+lci)=Tb(end/2+1:end);
end

% figure(201)
% imagesc(angle(coeff_matrix));
% figure(202)
% imagesc(angle(incident_wave_matrix));
% figure(203)
% imagesc(abs(coeff_matrix));
% figure(204)
% imagesc(abs(incident_wave_matrix));
