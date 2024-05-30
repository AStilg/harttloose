function [T,R,RgQ,Q]=tmatrix_cylinder_ebcm_romberg(nmax,k_medium,k_particle,wh);
% tmatrix_cylinder_ebcm - computes EBCM T-matrix for cylinders using
%                         romberg integration (it seems faster)
%
% usage:
% [T,R,RgQ,Q]=tmatrix_cylinder_ebcm_romberg(nmax,k_medium,k_particle,wh);
%
% where,
% nmax       - the number of "radial" basis functions to terminate. Found using
%              Wiscombe or Brock condition.
% k_medium   - wave number in the medium.
% k_particle - wave number in the particle.
% wh         - [width,height]
%
% T         - T-matrix
% R         - external reg-to-internal reg matrix.
%
% this code is needed to generate example outputs.
% 
% PACKAGE INFO

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
match_points=2*(2*nmax+2)+10000;
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

% oneeta=sqrt(1-eta.^2);
for ii=1:length(uniq_n)
    find_ni=find(n==ii);
    for jj=1:length(uniq_n)
        find_nj=find(n==jj);
        for kk=0:min(ii,jj)
            i_indx=find(m(find_ni)==kk);
            j_indx=find(m(find_nj)==kk);

            J11(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            J22(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            J12(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            J21(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            
            RgJ11(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            RgJ22(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            RgJ12(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            RgJ21(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),rtp(:,2),3);
        end
    end
end

for ii=1:length(uniq_n)
    find_ni=find(n==ii);
    for jj=1:length(uniq_n)
        find_nj=find(n==jj);
        for kk=1:min(ii,jj)
            i_indx=find(m(find_ni)==-kk);
            j_indx=find(m(find_nj)==-kk);

            mJ11(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            mJ22(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            mJ12(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N1(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            mJ21(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M1(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            
            mRgJ11(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(M3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            mRgJ22(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(N3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            mRgJ12(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(M3(find_ni(i_indx),:).',[],3),reshape(N(find_nj(j_indx),:)',[],3),rtp(:,2),3);
            mRgJ21(jj,ii,kk+1)=ndotxietaphicross_romberg(nxrtp,reshape(N3(find_ni(i_indx),:).',[],3),reshape(M(find_nj(j_indx),:)',[],3),rtp(:,2),3);
        end
    end
end
%%
Q11=zeros((nmax*(nmax+2)));
Q21=Q11;
Q12=Q11;
Q22=Q11;
RgQ11=Q11;
RgQ21=Q11;
RgQ12=Q11;
RgQ22=Q11;

for m=0:size(J11,3)-1
    m0=max(m,1);
    ci=combined_index([m0:nmax]',m);
    [CI,CJ]=ndgrid(ci,ci);
    Q11_index_map=sub2ind((nmax*(nmax+2))*[1,1],CI(:),CJ(:));
    
    Q11(Q11_index_map)=k_medium*(k_particle*J21([m0:nmax],[m0:nmax],m+1)+k_medium*J12([m0:nmax],[m0:nmax],m+1));
    Q21(Q11_index_map)=k_medium*(k_particle*J22([m0:nmax],[m0:nmax],m+1)+k_medium*J11([m0:nmax],[m0:nmax],m+1));
    Q12(Q11_index_map)=k_medium*(k_particle*J11([m0:nmax],[m0:nmax],m+1)+k_medium*J22([m0:nmax],[m0:nmax],m+1));
    Q22(Q11_index_map)=k_medium*(k_particle*J12([m0:nmax],[m0:nmax],m+1)+k_medium*J21([m0:nmax],[m0:nmax],m+1));
    
    RgQ11(Q11_index_map)=k_medium*(k_particle*RgJ21([m0:nmax],[m0:nmax],m+1)+k_medium*RgJ12([m0:nmax],[m0:nmax],m+1));
    RgQ21(Q11_index_map)=k_medium*(k_particle*RgJ22([m0:nmax],[m0:nmax],m+1)+k_medium*RgJ11([m0:nmax],[m0:nmax],m+1));
    RgQ12(Q11_index_map)=k_medium*(k_particle*RgJ11([m0:nmax],[m0:nmax],m+1)+k_medium*RgJ22([m0:nmax],[m0:nmax],m+1));
    RgQ22(Q11_index_map)=k_medium*(k_particle*RgJ12([m0:nmax],[m0:nmax],m+1)+k_medium*RgJ21([m0:nmax],[m0:nmax],m+1));
    if m>0
        Q11_index_map2=sub2ind((nmax*(nmax+2))*[1,1],CI(:)-2*m,CJ(:)-2*m);
        Q11(Q11_index_map2)=k_medium*(k_particle*mJ21([m0:nmax],[m0:nmax],m+1)+k_medium*mJ12([m0:nmax],[m0:nmax],m+1));
        Q21(Q11_index_map2)=k_medium*(k_particle*mJ22([m0:nmax],[m0:nmax],m+1)+k_medium*mJ11([m0:nmax],[m0:nmax],m+1));
        Q12(Q11_index_map2)=k_medium*(k_particle*mJ11([m0:nmax],[m0:nmax],m+1)+k_medium*mJ22([m0:nmax],[m0:nmax],m+1));
        Q22(Q11_index_map2)=k_medium*(k_particle*mJ12([m0:nmax],[m0:nmax],m+1)+k_medium*mJ21([m0:nmax],[m0:nmax],m+1));
        
        RgQ11(Q11_index_map2)=k_medium*(k_particle*mRgJ21([m0:nmax],[m0:nmax],m+1)+k_medium*mRgJ12([m0:nmax],[m0:nmax],m+1));
        RgQ21(Q11_index_map2)=k_medium*(k_particle*mRgJ22([m0:nmax],[m0:nmax],m+1)+k_medium*mRgJ11([m0:nmax],[m0:nmax],m+1));
        RgQ12(Q11_index_map2)=k_medium*(k_particle*mRgJ11([m0:nmax],[m0:nmax],m+1)+k_medium*mRgJ22([m0:nmax],[m0:nmax],m+1));
        RgQ22(Q11_index_map2)=k_medium*(k_particle*mRgJ12([m0:nmax],[m0:nmax],m+1)+k_medium*mRgJ21([m0:nmax],[m0:nmax],m+1));
    end
end
%% create a code to block up EBCM code to subspace.

T=zeros(2*size(Q11));
R=T;

Q=conj(1i*[Q11,-Q12;-Q21,Q22]/2/pi);
RgQ=conj(1i*[RgQ11,-RgQ12;-RgQ21,RgQ22]/2/pi);

MC=mode_couplings(nmax,[inf,1]);

%stabilize
Q=MC.*Q;
RgQ=MC.*RgQ;

eyeQ=eye(size(Q));

for ii=1:size(Q11,2)*2
    jvec=find(MC(:,ii));
    R(jvec,ii)=-(Q(:,jvec)\eyeQ(:,ii));
end

% T=-(Q'\RgQ')';

Qt=Q.';
RgQt=RgQ.';

for ii=1:size(Q11,2)*2

    % ivec=find(MC(ii,:));
    % Tu(ii,ivec)=-RgQu(ii,ivec)/Qu(ivec,ivec);


    ivec=find(MC(ii,:));
    T(ii,ivec)=-(Qt(:,ivec)\RgQt(:,ii)).';

end

