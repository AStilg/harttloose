function [T,T2]=tmatrix_cylinder_pm(nmax,k_medium,k_particle,wh);
% tmatrix_cylinder_pm - computes PM T-matrix for cylinders.
%
% usage:
% [T,T2]=tmatrix_cylinder_pm(nmax,k_medium,k_particle,wh)
%
% where,
% nmax       - the number of "radial" basis functions to terminate. Found using
%              Wiscombe or Brock condition.
% k_medium   - wave number in the medium.
% k_particle - wave number in the particle.
% wh         - [width,height]
%
% T         - T-matrix
% T2         - external reg-to-internal reg matrix.
%
% this code is needed to generate example outputs.
% 
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

ac=wh/2;

% nmax=10;
% ac=[1,1.001]/2;
% k_medium=2*pi;
% k_particle=2*pi*1.2;

nmaxs=nmax;
if numel(nmaxs)==1
    nmaxs=[nmax,nmax];
end
nmax=max(nmaxs);

syms theta phi real;
syms r real positive;
rv=r*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];

r1=subs(rv,r,ac(2)/cos(theta));
der1=diff(r1,theta);
dpr1=diff(r1,phi);
cp1=simplify(cross(der1,dpr1),'Steps',100);
% subs(subs(cp1,eta,.7),phi,0)

r2=subs(rv,r,ac(1)/sin(theta)/cos(phi));
der2=diff(r2,theta);
dpr2=diff(r2,phi);
cp2=simplify(cross(der2,dpr2));

r3=subs(rv,r,-ac(2)/cos(theta));
der3=diff(r3,theta);
dpr3=diff(r3,phi);
cp3=simplify(cross(der3,dpr3),'Steps',100);

% thetasol=isolate(r3(1)==r2(1),theta);


thetax = mod(atan2(ac(1),ac(2))+2*pi,2*pi);
% thetax-double(vpa(simplify(rhs(subs(thetasol,[phi],[0]))),16))
% c=1e-5;
ranges=[0,thetax];
ranges=[ranges;[ranges(2),pi-ranges(2)];[pi-ranges(2),pi]];

% create match points
match_points=2*(2*nmax+2)+1000;
match_points_pow2=round(log(match_points)/log(2));
match_points=2^(match_points_pow2)+1;

etan=zeros(match_points,size(ranges,1));

for ii=1:size(etan,2)
    etan(:,ii)=linspace(ranges(ii,1),ranges(ii,2),match_points)';
end

%% 1
fxyz1=matlabFunction(r1(1:2),'vars',[theta,phi]);
xyz1=fxyz1((etan(:,1)),0*etan(:,1));
xyz1(:,3)=double(r1(3));

fcp1=matlabFunction(cp1(3),'vars',[theta,phi]);
nxyz1=[0,0,1].*fcp1((etan(:,1)),0*etan(:,1));
%% 2
fxyz2=matlabFunction(r2(2:3),'vars',[theta,phi]); %simplify didn't work...
xyz2=fxyz2((etan(:,2)),0*etan(:,2));
xyz2=[double(r2(1))*ones(size(xyz2,1),1),xyz2];

fcp2=matlabFunction(cp2(1),'vars',[theta,phi]);
nxyz2=[1,0,0].*fcp2((etan(:,2)),0*etan(:,2));
%% 3
fxyz3=matlabFunction(r3(1:2),'vars',[theta,phi]); %simplify didn't work...
xyz3=fxyz3((etan(:,3)),0*etan(:,3));
xyz3(:,3)=double(r3(3));

fcp3=matlabFunction(cp3(3),'vars',[theta,phi]);
nxyz3=[0,0,1].*fcp3((etan(:,3)),0*etan(:,3));
%%


[nxrtp,rtp]=xyzv2rtpv([nxyz1;nxyz2;nxyz3],[xyz1;xyz2;xyz3]);
etan=etan(:);

% % plotting to check normals.
% figure(4)
% plot(xyz1(:,1),xyz1(:,3));
% hold on
% plot(xyz2(:,1),xyz2(:,3));
% plot(xyz3(:,1),xyz3(:,3));
% 
% plot((xyz1(:,1)+nxyz1(:,1).*[0,1])',(xyz1(:,3)+nxyz1(:,3).*[0,1])');
% hold on
% plot((xyz2(:,1)+nxyz2(:,1).*[0,1])',(xyz2(:,3)+nxyz2(:,3).*[0,1])');
% plot((xyz3(:,1)+nxyz3(:,1).*[0,1])',(xyz3(:,3)+nxyz3(:,3).*[0,1])');
% hold off
% %%
% 
% figure(5)
% plot(etan(:),rtp)
% 
% pause

%% all coordinates have been determined.

% we need evey element:
[n,m]=combined_index([1:nmax*(nmax+2)]');

[M1,N1,~,~,M,N]=vswf_vector(n,m,k_medium*rtp(:,1),rtp(:,2),rtp(:,3));
[M3,N3,~,~,~,~]=vswf_vector(n,m,k_particle*rtp(:,1),rtp(:,2),rtp(:,3),3);
M1=M1.';
N1=N1.';
M=M.';
N=N.';
M3=M3.';
N3=N3.';
for ii=1:size(M1,2)
    M1(:,ii)=reshape(cross(nxrtp,reshape(M1(:,ii),[],3)),[],1);
    N1(:,ii)=reshape(cross(nxrtp,reshape(N1(:,ii),[],3)),[],1);
    M(:,ii)=reshape(cross(nxrtp,reshape(M(:,ii),[],3)),[],1);
    N(:,ii)=reshape(cross(nxrtp,reshape(N(:,ii),[],3)),[],1);
end

for ii=1:size(M3,2)
    M3(:,ii)=reshape(cross(nxrtp,reshape(M3(:,ii),[],3)),[],1);
    N3(:,ii)=reshape(cross(nxrtp,reshape(N3(:,ii),[],3)),[],1);
end

CM1=[[M;N],[N;M]];
CM2=-[[M1;N1],[N1;M1]];
CM3=[[M3;k_particle/k_medium*N3],[N3;k_particle/k_medium*M3]];

symmetry=[inf,1];
modeMatrix=mode_couplings(nmaxs(1),symmetry);
modeMatrixI=mode_couplings(nmaxs(2),symmetry);

lci=size(modeMatrixI,2)/2;

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

