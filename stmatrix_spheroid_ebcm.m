function [T,R,c,RgQ,Q]=stmatrix_spheroid_ebcm(nmax,k_medium,k_particle,ac);
% stmatrix_spheroid_ebcm - computes EBCM sT-matrix for spheroids.
%
% usage:
% [sT,sT2,c_d,sRgQ,sQ,sTu]=stmatrix_spheroid_ebcm(nmax,k_medium,k_particle,ac);
%
% where,
% nmax       - the number of "radial" basis functions to terminate. Found using
%              Wiscombe or Brock condition at a maximum. Likely nmax will be less.
% k_medium   - wave number in the medium.
% k_particle - wave number in the particle.
% ac         - [semi-major,semi-minor] distance (oblate), swap for prolate.
%
% outputs,
% sT    - Tmatrix.
% sT2   - RgQ matrix.
% c_d   - interfocal distance (in units of ac).
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

% r=.5;
% ar=10;
% ii=1;
% ac=[ar(ii).^(-1/3),ar(ii).^(2/3)]*r;
% 
% nmax=ka2nmax(2*pi*1.2*r)+2;
% k_medium=2*pi;
% k_particle=2*pi*1.2;
% % ac=[.5/10,.5];

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

[eta,w]=gausslegendreroot(match_points);
eta=eta(:);

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
%% all coordinates have been determined.

% we need evey element:
[n,m]=combined_index([0:nmax*(nmax+2)]');

xep=[xi.*ones(size(eta)),eta,zeros(size(eta))];

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

oneeta=ones(size(eta));

h=spheroidal_scale_factors(isProlate,k_medium*c,xep(:,1),xep(:,2));

dS=w.*h(:,2).*h(:,3); %include scale factors and GL weights

% oneeta=sqrt(1-eta.^2);
for ii=1:length(uniq_n)
    find_ni=find(n==ii-1);
    for jj=1:length(uniq_n)
        find_nj=find(n==jj-1);
        for kk=0:min(ii-1,jj-1)
            i_indx=find(m(find_ni)==kk);
            j_indx=find(m(find_nj)==kk);

            J11(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),dS);
            J22(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),dS);
            J12(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),dS);
            J21(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),dS);
            
            RgJ11(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),dS);
            RgJ22(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),dS);
            RgJ12(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),dS);
            RgJ21(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),dS);
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

            mJ11(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),dS);
            mJ22(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),dS);
            mJ12(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),dS);
            mJ21(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),dS);
            
            mRgJ11(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),dS);
            mRgJ22(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),dS);
            mRgJ12(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),dS);
            mRgJ21(jj,ii,kk+1)=ndotxietaphicross([1,0,0].*oneeta,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),dS);
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
k_m2=(k_medium);
k_p2=(k_particle);
% k_medium=k_particle;
% k_particle=k_m2;
for m=0:size(J11,3)-1
    ci=combined_index([m:nmax]',m)+1;
    [CI,CJ]=ndgrid(ci,ci);
    Q11_index_map=sub2ind((nmax*(nmax+2)+1)*[1,1],CI(:),CJ(:));
    
    Q11(Q11_index_map)=k_m2*(k_p2*J21([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*J12([m:nmax]+1,[m:nmax]+1,m+1));
    Q21(Q11_index_map)=k_m2*(k_p2*J22([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*J11([m:nmax]+1,[m:nmax]+1,m+1));
    Q12(Q11_index_map)=k_m2*(k_p2*J11([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*J22([m:nmax]+1,[m:nmax]+1,m+1));
    Q22(Q11_index_map)=k_m2*(k_p2*J12([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*J21([m:nmax]+1,[m:nmax]+1,m+1));
    
    RgQ11(Q11_index_map)=k_m2*(k_p2*RgJ21([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*RgJ12([m:nmax]+1,[m:nmax]+1,m+1));
    RgQ21(Q11_index_map)=k_m2*(k_p2*RgJ22([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*RgJ11([m:nmax]+1,[m:nmax]+1,m+1));
    RgQ12(Q11_index_map)=k_m2*(k_p2*RgJ11([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*RgJ22([m:nmax]+1,[m:nmax]+1,m+1));
    RgQ22(Q11_index_map)=k_m2*(k_p2*RgJ12([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*RgJ21([m:nmax]+1,[m:nmax]+1,m+1));
    if m>0
        Q11_index_map2=sub2ind((nmax*(nmax+2)+1)*[1,1],CI(:)-2*m,CJ(:)-2*m);
        Q11(Q11_index_map2)=k_m2*(k_p2*mJ21([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*mJ12([m:nmax]+1,[m:nmax]+1,m+1));
        Q21(Q11_index_map2)=k_m2*(k_p2*mJ22([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*mJ11([m:nmax]+1,[m:nmax]+1,m+1));
        Q12(Q11_index_map2)=k_m2*(k_p2*mJ11([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*mJ22([m:nmax]+1,[m:nmax]+1,m+1));
        Q22(Q11_index_map2)=k_m2*(k_p2*mJ12([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*mJ21([m:nmax]+1,[m:nmax]+1,m+1));
        
        RgQ11(Q11_index_map2)=k_m2*(k_p2*mRgJ21([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*mRgJ12([m:nmax]+1,[m:nmax]+1,m+1));
        RgQ21(Q11_index_map2)=k_m2*(k_p2*mRgJ22([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*mRgJ11([m:nmax]+1,[m:nmax]+1,m+1));
        RgQ12(Q11_index_map2)=k_m2*(k_p2*mRgJ11([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*mRgJ22([m:nmax]+1,[m:nmax]+1,m+1));
        RgQ22(Q11_index_map2)=k_m2*(k_p2*mRgJ12([m:nmax]+1,[m:nmax]+1,m+1)+k_m2*mRgJ21([m:nmax]+1,[m:nmax]+1,m+1));
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
    % Tu(ivec,ii)=-Qu(ivec,ivec)\RgQu(ivec,ii);


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
