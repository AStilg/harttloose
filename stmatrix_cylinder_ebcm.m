function [T,R,c,RgQ,Q]=stmatrix_cylinder_ebcm(nmax,k_medium,k_particle,wh);
% stmatrix_cylinder_ebcm - computes EBCM sT-matrix for cylinders.
%
% usage:
% [sT,sT2,c_d,sRgQ,sQ]=stmatrix_cylinder_ebcm(nmax,k_medium,k_particle,wh);
%
% where,
% nmax       - the number of "radial" basis functions to terminate.
% k_medium   - wave number in the medium.
% k_particle - wave number in the particle.
% wh         - [width,height]
%
% outputs,
% sT    - Tmatrix.
% sT2   - RgQ matrix.
% c_d   - interfocal distance (in units of wh).
% sQ    - regular internal-to-regular external couplings
% sRgQ  - regular internal-to-outgoing external couplings
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
%etax-double(vpa(rhs(subs(etasol,[phi,c_sym],[0,c])),16))
% c=1e-5;
ranges=[-1,-double(vpa(rhs(subs(etasol,[phi,c_sym],[0,c])),16))];
ranges=[ranges;[ranges(2),-ranges(2)];[-ranges(2),1]];

% create match points
match_points=2*(2*nmax+2)+100;

etan=zeros(match_points,size(ranges,1));
w=etan;

for ii=1:size(etan,2)
    [etan(:,ii),w(:,ii)]=gausslegendreroot(match_points,ranges(ii,:));
end

%% 1
fxyz1=matlabFunction(r1(1:2),'vars',[c_sym,eta,phi]);
xyz1=fxyz1(c,etan(:,1),0*etan(:,1));
xyz1(:,3)=double(r1(3));

fcp1=matlabFunction(cp1(3),'vars',[c_sym,eta,phi]);
nxyz1=[0,0,1].*fcp1(c,etan(:,1),0*etan(:,1));
%% 2
fxyz2=matlabFunction(r2,'vars',[c_sym,eta,phi]); %simplify didn't work...
xyz2=fxyz2(c,etan(:,2),0*etan(:,2));

fcp2=matlabFunction(cp2(1),'vars',[c_sym,eta,phi]);
nxyz2=[1,0,0].*fcp2(c,etan(:,2),0*etan(:,2));
%% 3
fxyz3=matlabFunction(r3(1:2),'vars',[c_sym,eta,phi]); %simplify didn't work...
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

uniq_n=unique(n);
%positive m;
J11=zeros(length(uniq_n),length(uniq_n),length(max(uniq_n))+1);
J22=J11;
J12=J11;
J21=J11;
RgJ11=J11;
RgJ22=J11;
RgJ12=J11;
RgJ21=J11;
mJ11=J11;
mJ12=J11;
mJ21=J11;
mJ22=J11;
mRgJ11=J11;
mRgJ12=J11;
mRgJ21=J11;
mRgJ22=J11;

dS=w;

% oneeta=sqrt(1-eta.^2);
for ii=1:length(uniq_n)
    find_ni=find(n==ii-1);
    for jj=1:length(uniq_n)
        find_nj=find(n==jj-1);
        for kk=0:min(ii-1,jj-1)
            i_indx=find(m(find_ni)==kk);
            j_indx=find(m(find_nj)==kk);

            J11(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),dS);
            J22(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),dS);
            J12(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),dS);
            J21(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),dS);
            
            RgJ11(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),dS);
            RgJ22(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),dS);
            RgJ12(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),dS);
            RgJ21(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),dS);
        end
    end
end

for ii=1:length(uniq_n)
    find_ni=find(n==ii-1);
    for jj=1:length(uniq_n)
        find_nj=find(n==jj-1);
        for kk=0:min(ii-1,jj-1)
            i_indx=find(m(find_ni)==-kk);
            j_indx=find(m(find_nj)==-kk);

            mJ11(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),dS);
            mJ22(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),dS);
            mJ12(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),dS);
            mJ21(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),dS);
            
            mRgJ11(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),dS);
            mRgJ22(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),dS);
            mRgJ12(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),dS);
            mRgJ21(jj,ii,kk+1)=ndotxietaphicross(nxepv,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),dS);
        end
    end
end
%%
Q11=zeros((nmax*(nmax+2)+1));
Q21=Q11;
Q12=Q11;
Q22=Q11;
RgQ11=Q11;
RgQ21=Q11;
RgQ12=Q11;
RgQ22=Q11;
J11m=Q11;
J12m=Q11;
% k_m2=k_medium;
% k_medium=k_particle;
% k_particle=k_m2;
for m=0:size(J11,3)-1
    ci=combined_index([m:nmax]',m)+1;
    [CI,CJ]=ndgrid(ci,ci);
    Q11_index_map=sub2ind((nmax*(nmax+2)+1)*[1,1],CI(:),CJ(:));
    
    Q11(Q11_index_map)=k_medium*(k_particle*J21([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*J12([m:nmax]+1,[m:nmax]+1,m+1));
    Q21(Q11_index_map)=k_medium*(k_particle*J22([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*J11([m:nmax]+1,[m:nmax]+1,m+1));
    Q12(Q11_index_map)=k_medium*(k_particle*J11([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*J22([m:nmax]+1,[m:nmax]+1,m+1));
    Q22(Q11_index_map)=k_medium*(k_particle*J12([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*J21([m:nmax]+1,[m:nmax]+1,m+1));
    
    RgQ11(Q11_index_map)=k_medium*(k_particle*RgJ21([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*RgJ12([m:nmax]+1,[m:nmax]+1,m+1));
    RgQ21(Q11_index_map)=k_medium*(k_particle*RgJ22([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*RgJ11([m:nmax]+1,[m:nmax]+1,m+1));
    RgQ12(Q11_index_map)=k_medium*(k_particle*RgJ11([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*RgJ22([m:nmax]+1,[m:nmax]+1,m+1));
    RgQ22(Q11_index_map)=k_medium*(k_particle*RgJ12([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*RgJ21([m:nmax]+1,[m:nmax]+1,m+1));
    if m>0
        Q11_index_map2=sub2ind((nmax*(nmax+2)+1)*[1,1],CI(:)-2*m,CJ(:)-2*m);
        Q11(Q11_index_map2)=k_medium*(k_particle*mJ21([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*mJ12([m:nmax]+1,[m:nmax]+1,m+1));
        Q21(Q11_index_map2)=k_medium*(k_particle*mJ22([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*mJ11([m:nmax]+1,[m:nmax]+1,m+1));
        Q12(Q11_index_map2)=k_medium*(k_particle*mJ11([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*mJ22([m:nmax]+1,[m:nmax]+1,m+1));
        Q22(Q11_index_map2)=k_medium*(k_particle*mJ12([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*mJ21([m:nmax]+1,[m:nmax]+1,m+1));
        
        RgQ11(Q11_index_map2)=k_medium*(k_particle*mRgJ21([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*mRgJ12([m:nmax]+1,[m:nmax]+1,m+1));
        RgQ21(Q11_index_map2)=k_medium*(k_particle*mRgJ22([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*mRgJ11([m:nmax]+1,[m:nmax]+1,m+1));
        RgQ12(Q11_index_map2)=k_medium*(k_particle*mRgJ11([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*mRgJ22([m:nmax]+1,[m:nmax]+1,m+1));
        RgQ22(Q11_index_map2)=k_medium*(k_particle*mRgJ12([m:nmax]+1,[m:nmax]+1,m+1)+k_medium*mRgJ21([m:nmax]+1,[m:nmax]+1,m+1));
    end
end
%% create a code to block up EBCM code to subspace.

T=zeros(2*size(Q11));
Tu=T;
R=T;

Qu=conj(1i*[Q11,-Q12;-Q21,Q22]/2/pi);
RgQu=conj(1i*[RgQ11,-RgQ12;-RgQ21,RgQ22]/2/pi);
% eyeR=eye(size(RgQu));
MC=mode_couplings_spheroidal(nmax,[inf,1]);


mode_index=[1:nmax*(nmax+2)];
n=floor(sqrt(mode_index));
N2=sparse(mode_index,mode_index,n.*(n+1));

UNUd=spherical_to_spheroidal(isProlate,k_medium*c,N2);
UUNNUU=blkdiag(UNUd,UNUd);
UiNUd=spherical_to_spheroidal(isProlate,k_medium*c,inv(N2));
UUiNNUU=blkdiag(UiNUd,UiNUd);

Q=UUiNNUU*Qu;
RgQ=UUiNNUU*RgQu;

for ii=1:size(Q11,2)*2
    jvec=find(MC(:,ii));
    R(jvec,ii)=-(Qu(:,jvec)\UUNNUU(:,ii));
end

Qut=Qu.';
RgQut=RgQu.';


for ii=1:size(Q11,2)*2

    % ivec=find(MC(ii,:));
    % Tu(ii,ivec)=-RgQu(ii,ivec)/Qu(ivec,ivec);


    ivec=find(MC(ii,:));
    Tu(ii,ivec)=-(Qut(:,ivec)\RgQut(:,ii)).';

end


T=UUiNNUU*Tu*UUNNUU;

% for ii=1:size(Q11,2)*2
%     jvec=find(MC(:,ii));
%     T(jvec,ii)=RgQ(jvec,jvec)*R(jvec,ii)/2; %transpose
% end

% % Trial code"
% for ii=1:size(Q11,2)*2
%     jvec=find(MC(ii,:));
%     T(jvec,ii)=-Q(:,jvec)\RgQ(:,ii); %transpose
%     R(jvec,ii)=-2*(eyeR(:,jvec)\RgQ(:,ii));
% end
