function [T,T2,coeff_matrix,incident_wave_matrix]=tmatrix_spheroid_pm(nmax,k_medium,k_particle,ac);
% tmatrix_spheroid_pm - computes point-matched T-matrix for spheroids
% using point-matching method.
%
% usage:
% [T,T2]=tmatrix_spheroid_pm(nmax,k_medium,k_particle,ac);
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
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

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

% create match points
match_points=2*(2*nmax+2)+100;

% eta=linspace(-1+1/match_points,1-1/match_points,match_points);
eta=gausslegendreroot(match_points);
% eta=cumsum(rand(1,match_points)*4/match_points)-1;
eta=eta(:);

eta(eta<0)=[];
match_points=length(eta);

% eta is costheta
xi=sqrt(1-eta.^2); %not xi
xyz=[xi*minor_axis,0*xi,eta*major_axis];

nxyz=[aspect_ratio*xi*minor_axis,0*eta,eta*major_axis/aspect_ratio]./sqrt((eta/aspect_ratio*major_axis).^2+(aspect_ratio*xi*minor_axis).^2);

xy = sqrt( xyz(:,1).^2 + xyz(:,2).^2 );
theta = mod(atan2(xy,xyz(:,3))+2*pi,2*pi);
phi = mod(atan2(xyz(:,2),xyz(:,1))+2*pi,2*pi);
r = sqrt( xy.^2 + xyz(:,3).^2 );
rtp=[r(:),theta(:),phi(:)];

J=[sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta);...
    cos(theta).*cos(phi),cos(theta).*sin(phi),-sin(theta);...
    -sin(phi),cos(phi),zeros(size(theta))];

%Separate the Jacobian and multiply for each unit vector.
normals = [dot(J(1:length(theta),:),nxyz,2),dot(J(length(theta)+1:2*length(theta),:),nxyz,2),dot(J(2*length(theta)+1:3*length(theta),:),nxyz,2)];

nr=normals(:,1);
nt=normals(:,2);
np=normals(:,3);
r=rtp(:,1);
t=rtp(:,2);
p=rtp(:,3);
%% pack matrix of fields:

total_orders=nmax*(nmax+2);
ci_total=[1:total_orders]';
lci=length(ci_total);

%stick with E_perp continuous opposed to D_para continuous. return to
%check.
% coeff_matrix = zeros(6*match_points,4*lci);
% incident_wave_matrix = zeros(6*match_points,2*lci);

[n,m]=combined_index(ci_total);

[M2t,N2t,~,~,M1t,N1t]=vswf_vector(n,m,k_medium*r,t,p);
[M3t,N3t]=vswf_vector(n,m,k_particle*r,t,p,3);
M2t=M2t.';
N2t=N2t.';
M1t=M1t.';
N1t=N1t.';
M3t=M3t.';
N3t=N3t.';
for ii=1:size(M1t,2)
    M1t(:,ii)=reshape(cross([nr,nt,np],reshape(M1t(:,ii),[],3)),[],1);
    N1t(:,ii)=reshape(cross([nr,nt,np],reshape(N1t(:,ii),[],3)),[],1);
    M2t(:,ii)=reshape(cross([nr,nt,np],reshape(M2t(:,ii),[],3)),[],1);
    N2t(:,ii)=reshape(cross([nr,nt,np],reshape(N2t(:,ii),[],3)),[],1);
end

for ii=1:size(M3t,2)
    M3t(:,ii)=reshape(cross([nr,nt,np],reshape(M3t(:,ii),[],3)),[],1);
    N3t(:,ii)=reshape(cross([nr,nt,np],reshape(N3t(:,ii),[],3)),[],1);
end

CM1=[[M1t;N1t],[N1t;M1t]];
CM2=-[[M2t;N2t],[N2t;M2t]];
CM3=[[M3t;k_particle/k_medium*N3t],[N3t;k_particle/k_medium*M3t]];

symmetry=[inf,1];
modeMatrix=mode_couplings(nmaxs(1),symmetry);
modeMatrixI=mode_couplings(nmaxs(2),symmetry);


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
    
%     Ta=[CM2(:,index),CM3(:,index)]\(CM1(:,ii)); %E fields
%     Tb=[CM2(:,index2),CM3(:,index2)]\CM1(:,ii+lci); %H fields
    
    
    Ta=([CM2(:,index),CM3(:,index)]'*[CM2(:,index),CM3(:,index)])\([CM2(:,index),CM3(:,index)]'*(CM1(:,ii))); %E fields
    Tb=([CM2(:,index2),CM3(:,index2)]'*[CM2(:,index2),CM3(:,index2)])\([CM2(:,index2),CM3(:,index2)]'*CM1(:,ii+lci)); %H fields
    
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
