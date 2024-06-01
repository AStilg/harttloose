function [T,T2,c]=stmatrix_cylinder_pm(nmax,k_medium,k_particle,wh);
% stmatrix_cylinder_pm - computes sT-matrix for cylinders
%                       using point-matching method.
%
% usage:
% [T,T2,c_d]=stmatrix_cylinder_pm(nmax,k_medium,k_particle,wh);
%
% where,
% nmax       - the number of "radial" basis functions to terminate.
% k_medium   - wave number in the medium.
% k_particle - wave number in the particle.
% wh         - [width,height] distance (oblate).
%
% outputs,
% T    - Tmatrix.
% T2   - RgQ matrix.
% c_d  - interfocal distance (in units of wh).
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

ac=wh/2;

nmaxs=nmax;
if numel(nmaxs)==1
    nmaxs=[nmax,nmax];
end
nmax=max(nmaxs);

k_rel=k_particle/k_medium;

isProlate=ac(2)>ac(1);
sigma=2*isProlate-1;


syms eta phi real;
syms c_sym xi real positive;
r=c_sym*[sqrt(xi.^2-sigma).*sqrt(1-eta.^2).*cos(phi),sqrt(xi.^2-sigma).*sqrt(1-eta.^2).*sin(phi),xi.*eta];

r1=subs(r,xi,-ac(2)/eta/c_sym);
der1=diff(r1,eta);
dpr1=diff(r1,phi);
cp1=simplify(cross(der1,dpr1),'Steps',100);
% subs(subs(cp1,eta,.7),phi,0)

r2=subs(r,xi,sqrt(ac(1)^2/(1-eta^2)/c_sym^2/cos(phi)^2+sigma));
der2=diff(r2,eta);
dpr2=diff(r2,phi);
cp2=simplify(cross(der2,dpr2));

r3=subs(r,xi,ac(2)/eta/c_sym);
der3=diff(r3,eta);
dpr3=diff(r3,phi);
cp3=simplify(cross(der3,dpr3),'Steps',100);

etasol=isolate(r3(1)==ac(1),eta);
c=(sqrt(sigma*((ac(2)/ac(1))^2-1))*ac(1));


[xix,etax,phix]=xyz2xietaphi(isProlate,c,[ac(1),0,ac(2)]);
% etax-double(vpa(rhs(subs(etasol,[phi,c_sym],[0,c])),16));
% c=1e-5;
ranges=[-1,-double(vpa(rhs(subs(etasol,[phi,c_sym],[0,c])),16))];
ranges=[ranges;[ranges(2),-ranges(2)];[-ranges(2),1]];

% create match points
match_points=2*(2*nmax+2)+10;

etan=zeros(match_points,size(ranges,1));
w=etan;

for ii=1:size(etan,2)
    [etan(:,ii),w(:,ii)]=gausslegendreroot(match_points,ranges(ii,:));
end

%% 1
fxyz1=matlabFunction(r1(1:2));
xyz1=fxyz1(c,etan(:,1),0*etan(:,1));
xyz1(:,3)=double(r1(3));

fcp1=matlabFunction(cp1(3),'vars',[c_sym,eta,phi]);
nxyz1=[0,0,1].*fcp1(c,etan(:,1),0*etan(:,1));
%% 2
fxyz2=matlabFunction(r2); %simplify didn't work...
xyz2=fxyz2(c,etan(:,2),0*etan(:,2));

fcp2=matlabFunction(cp2(1),'vars',[c_sym,eta,phi]);
nxyz2=[1,0,0].*fcp2(c,etan(:,2),0*etan(:,2));
%% 3
fxyz3=matlabFunction(r3(1:2)); %simplify didn't work...
xyz3=fxyz3(c,etan(:,3),0*etan(:,3));
xyz3(:,3)=double(r3(3));

fcp3=matlabFunction(cp3(3),'vars',[c_sym,eta,phi]);
nxyz3=[0,0,1].*fcp3(c,etan(:,3),0*etan(:,3));
%%

etan=etan(:);
w=w(:);

[nxepv,xep]=xyzv2xietaphiv(isProlate,c,-[nxyz1;nxyz2;nxyz3],[xyz1;xyz2;xyz3]);

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
% plot(etan(:),xep)
% 
% pause

%% all coordinates have been determined.

% we need evey element:
[n,m]=combined_index([0:nmax*(nmax+2)]');

[M1,N1,~,~,M,N]=spheroidalvwf(isProlate,n,m,k_medium*c,xep(:,1:3));
[M3,N3,~,~,~,~]=spheroidalvwf(isProlate,n,m,k_particle*c,xep(:,1:3),3);

M1=M1.';
N1=N1.';
M=M.';
N=N.';
M3=M3.';
N3=N3.';

for ii=1:size(M1,2)
    M1(:,ii)=reshape(cross(nxepv,reshape(M1(:,ii),[],3)),[],1);
    N1(:,ii)=reshape(cross(nxepv,reshape(N1(:,ii),[],3)),[],1);
    M(:,ii)=reshape(cross(nxepv,reshape(M(:,ii),[],3)),[],1);
    N(:,ii)=reshape(cross(nxepv,reshape(N(:,ii),[],3)),[],1);
end

for ii=1:size(M3,2)
    M3(:,ii)=reshape(cross(nxepv,reshape(M3(:,ii),[],3)),[],1);
    N3(:,ii)=reshape(cross(nxepv,reshape(N3(:,ii),[],3)),[],1);
end

CM1=[[M;N],[N;M]];
CM2=-[[M1;N1],[N1;M1]];
CM3=[[M3;k_particle/k_medium*N3],[N3;k_particle/k_medium*M3]];

symmetry=[inf,1];
modeMatrix=mode_couplings_spheroidal(nmaxs(1),symmetry);
modeMatrixI=mode_couplings_spheroidal(nmaxs(2),symmetry);

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

