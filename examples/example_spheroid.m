% example_spheroid.m
%
% This code was successfully run on a computer with: Intel Core i7 7700K, 32GB ram in MATLAB R2020b.
%
% 

addpath ../ ../provided/ ./helpers/

%all distance units are \lambda.
k_medium=2*pi;
k_particle=2*pi*1.2;

ac=[1.0,4]/.5; %8\lambda semi-major axis, 2\lambda semi-minor axis.

% set nmax
nmaxmin=ka2nmax(2*pi*ac*k_particle/k_medium);
nmax=ceil(mean(nmaxmin));

% calculate sT-matrix with pm and ebcm.
warning('off','MATLAB:rankDeficientMatrix'); %cull warnings... there are many.
isProlate=true;
tic
[sT_pm,sR_pm,c]=stmatrix_spheroid_pm(nmax,k_medium,k_particle,ac);
t_stmatrix_pm=toc
tic
[T_pm2,R_pm2]=tmatrix_spheroid_pm(nmax,k_medium,k_particle,ac); % this file is a particular implementation for this toolbox. It appears nowhere else.
t_tmatrix_pm=toc
tic
[sT_ebcm,sR_ebcm,c]=stmatrix_spheroid_ebcm(nmax,k_medium,k_particle,ac);
t_stmatrix_ebcm=toc
warning('on','MATLAB:rankDeficientMatrix');

%% NOTE: Code in this section requires ott from github and smarties 1.01
% %smarties T-matrix
% T_sm=ott.TmatrixSmarties(ac(1),ac(2),'wavelength0',1,'k_medium',k_medium,'k_particle',k_particle,'nmax',nmaxmin(2));;
% %beam
% ab=ott.BscPmGauss('lg',[0,0],'angle',64.5,'polarisation',[1,1i],'angular_scaling','sintheta');
%T_sm=T_sm.data;
%ab=[ab.a;ab.b];
%% Instead we load some we made earlier:
load('example_spheroid_data.mat');
%%
%create spheroidal beam coefficients
sa=spherical_to_spheroidal(isProlate,k_medium*c,ab(1:end/2),nmax);
sb=spherical_to_spheroidal(isProlate,k_medium*c,ab(end/2+1:end),nmax);

%create spheroidal sT_matrix and sR_matrix from spherical PM:
sT_pm2=[spherical_to_spheroidal(isProlate,k_medium*c,T_pm2(1:end/2,1:end/2),nmax),spherical_to_spheroidal(isProlate,k_medium*c,T_pm2(1:end/2,end/2+1:end),nmax);
        spherical_to_spheroidal(isProlate,k_medium*c,T_pm2(end/2+1:end,1:end/2),nmax),spherical_to_spheroidal(isProlate,k_medium*c,T_pm2(end/2+1:end,end/2+1:end),nmax)];
%the external to internal matrix has two mediums of operation:
sR_pm2=[spherical_to_spheroidal(isProlate,[k_particle*c,k_medium*c],R_pm2(1:end/2,1:end/2),nmax),spherical_to_spheroidal(isProlate,[k_particle*c,k_medium*c],R_pm2(1:end/2,end/2+1:end),nmax);
        spherical_to_spheroidal(isProlate,[k_particle*c,k_medium*c],R_pm2(end/2+1:end,1:end/2),nmax),spherical_to_spheroidal(isProlate,[k_particle*c,k_medium*c],R_pm2(end/2+1:end,end/2+1:end),nmax)];

% compute all scattered modes, spheroidal coefficients all start with "s"
spq_pm=sT_pm*[sa;sb];
scd_pm=sR_pm*[sa;sb];
spq_pm2=sT_pm2*[sa;sb];
scd_pm2=sR_pm2*[sa;sb];

%% now generate fields for these modes:
% first set up coordinates for all plotting:
r_v=max(ac)*1.05;
nx=100;
nz=100;

xv=linspace(-r_v,r_v,nx);
zv=linspace(-r_v,r_v,nz);
[X,Z,Y]=meshgrid(xv,zv,0);

%create major, minor and ellipsoids to draw:
th=linspace(-pi,pi);

ct=cos(th);
st=sin(th);

majx=ac(2)*ct;
majy=ac(2)*st;
minx=ac(1)*ct;
miny=ac(1)*st;
elx=ac(1)*ct;
ely=ac(2)*st;
lel=1;

%create coordinates for both sets of calculations:
XYZ=[X(:),Y(:),Z(:)];

%% spheroidal wavefunction calculations.
xep=xyz2xietaphi(isProlate,c,XYZ);


smodes_used=find(abs(sa))-1;
[snn,smm]=combined_index(smodes_used);

warning('The most time is spent calculating the fields. This also takes the most memory. If you are memory limited, you can break up the below into several steps. I''ve done it this way because it uses less code.')

[sM_sca,sN_sca,~,~,sM_inc,sN_inc]=spheroidalvwf(isProlate,snn,smm,k_medium*c,xep(:,1),xep(:,2),xep(:,3));
[sM_int,sN_int]=spheroidalvwf(isProlate,snn,smm,k_particle*c,xep(:,1),xep(:,2),xep(:,3),3);
%% compute s-coefficients.

sE_inc=reshape([sM_inc.',sN_inc.']*[sa(smodes_used+1);sb(smodes_used+1)],[],3);
% spheroidal pm first
offset=length(spq_pm)/2;
sE_sca_pm=reshape([sM_sca.',sN_sca.']*spq_pm([smodes_used+1;offset+(smodes_used+1)]),[],3);
sE_int_pm=reshape([sM_int.',sN_int.']*scd_pm([smodes_used+1;offset+(smodes_used+1)]),[],3);

% converted pm second
sE_sca_pm2=reshape([sM_sca.',sN_sca.']*spq_pm2([smodes_used+1;offset+(smodes_used+1)]),[],3);
sE_int_pm2=reshape([sM_int.',sN_int.']*scd_pm2([smodes_used+1;offset+(smodes_used+1)]),[],3);

% To make the fields "look" continuous for a dielectric we can weight the
% normal field component for the internal field, D=\eps E. \eps\propto k^2
sE_int_pm(:,1)=(k_particle/k_medium)^2*sE_int_pm(:,1);
sE_int_pm2(:,1)=(k_particle/k_medium)^2*sE_int_pm2(:,1);

%% plot modified spheroidal fields

mask=sqrt((X/ac(1)).^2+(Y/ac(1)).^2+(Z/ac(2)).^2)>1;
sE_pm=mask(:).*(sE_inc+sE_sca_pm)+(1-mask(:)).*(sE_int_pm);
sI_pm=reshape(sum(abs(sE_pm).^2,2),size(X));

sE_pm2=mask(:).*(sE_inc+sE_sca_pm2)+(1-mask(:)).*(sE_int_pm2);
sI_pm2=reshape(sum(abs(sE_pm2).^2,2),size(X));

cmapR=hsv(256);
cmapR(end,:)=1;

figure(1)
set(gcf,'position',[800         400         564         543])
colormap(cmapR)
subplot(2,2,1)
imagesc(xv,zv,abs(sI_pm))
hold on
plot(majx,majy,'w--');
plot(minx,miny,'w--');
plot(elx,ely,'w--','linewidth',lel);
xlabel('x [\lambda]')
ylabel('z [\lambda]')
grid on
hold off
axis equal square
axis(r_v*[-1,1,-1,1])
cax=caxis;
cvz=linspace(cax(1),cax(2)/30,2);
cax=cvz;
caxis(cax);
subplot(2,2,3)
imagesc(xv,zv,abs(sI_pm2))
hold on
plot(majx,majy,'w--');
plot(minx,miny,'w--');
plot(elx,ely,'w--','linewidth',lel);
xlabel('x [\lambda]')
ylabel('z [\lambda]')
grid on
hold off
axis equal square
axis(r_v*[-1,1,-1,1])
caxis(cax)

cxx=[-1.5,1.5];
figure(23)
subplot(1,2,1)
[C,h]=contourf(X,Z,log10(abs(sI_pm2-sI_pm)./abs(sI_pm)*100),linspace((cxx(1)),(cxx(2)),4));
% clabel(C,h)
hold on
plot(majx,majy,'w--');
plot(minx,miny,'w--');
plot(elx,ely,'w--','linewidth',lel*2);
xlabel('x [\lambda]')
ylabel('z [\lambda]')
hold off
% caxis([0,.5])
axis equal tight
h=colorbar
% colormap(parula(8))
grid on
set(gca,'ydir','reverse')
h.Label.String='log_{10} |percent difference|'


sE_pm=mask(:).*(sE_sca_pm)+(1-mask(:)).*(sE_int_pm-sE_inc);
sI_pm=reshape(sum(abs(sE_pm).^2,2),size(X));

sE_pm2=mask(:).*(sE_sca_pm2)+(1-mask(:)).*(sE_int_pm2-sE_inc);
sI_pm2=reshape(sum(abs(sE_pm2).^2,2),size(X));


subplot(1,2,2)
[C,h]=contourf(X,Z,log10(abs(sI_pm2-sI_pm)./abs(sI_pm)*100),linspace((cxx(1)),(cxx(2)),4));
% clabel(C,h)
hold on
plot(majx,majy,'w--');
plot(minx,miny,'w--');
plot(elx,ely,'w--','linewidth',lel*2);
xlabel('x [\lambda]')
ylabel('z [\lambda]')
hold off
% caxis([0,.5])
axis equal tight
h=colorbar
% colormap(parula(8))
grid on
set(gca,'ydir','reverse')
h.Label.String='log_{10} |percent difference|'

%% now for the spherical wavefunction version.
% first find the spherical T, and R, matrix:

offset=length(ab)/2;
ab2=ab;
ab2=ab2([[1:nmax*(nmax+2)],offset+[1:nmax*(nmax+2)]]);


T_pm=[spheroidal_to_spherical(isProlate,k_medium*c,sT_pm(1:end/2,1:end/2),nmax),spheroidal_to_spherical(isProlate,k_medium*c,sT_pm(1:end/2,end/2+1:end),nmax);
        spheroidal_to_spherical(isProlate,k_medium*c,sT_pm(end/2+1:end,1:end/2),nmax),spheroidal_to_spherical(isProlate,k_medium*c,sT_pm(end/2+1:end,end/2+1:end),nmax)];
R_pm=[spheroidal_to_spherical(isProlate,[k_particle*c,k_medium*c],sR_pm(1:end/2,1:end/2),nmax),spheroidal_to_spherical(isProlate,[k_particle*c,k_medium*c],sR_pm(1:end/2,end/2+1:end),nmax);
        spheroidal_to_spherical(isProlate,[k_particle*c,k_medium*c],sR_pm(end/2+1:end,1:end/2),nmax),spheroidal_to_spherical(isProlate,[k_particle*c,k_medium*c],sR_pm(end/2+1:end,end/2+1:end),nmax)];

%calculate scattered/internal coefficients

pq_pm=T_pm*[ab2];
cd_pm=R_pm*[ab2];

pq_pm2=T_pm2*[ab2];
cd_pm2=R_pm2*[ab2];

%% spherical wavefunctions... much faster.
modes_used=find(abs(ab2(1:end/2)));
[nn,mm]=combined_index(modes_used);

[r,t,p]=xietaphi2rtp(isProlate,c,xep(:,1),xep(:,2),xep(:,3));

[M_sca,N_sca,~,~,M_inc,N_inc]=vswf_vector(nn,mm,k_medium*r,t,p);
[M_int,N_int]=vswf_vector(nn,mm,k_particle*r,t,p,3);
%% compute fields in spherical wavefunctions.
offset=length(ab2)/2;

E_inc=reshape([M_inc.',N_inc.']*[ab2(modes_used);ab2(offset+modes_used)],[],3);
% spheroidal pm first

E_sca_pm=reshape([M_sca.',N_sca.']*pq_pm([modes_used;offset+(modes_used)]),[],3);
E_int_pm=reshape([M_int.',N_int.']*cd_pm([modes_used;offset+(modes_used)]),[],3);

% converted pm second
E_sca_pm2=reshape([M_sca.',N_sca.']*pq_pm2([modes_used;offset+(modes_used)]),[],3);
E_int_pm2=reshape([M_int.',N_int.']*cd_pm2([modes_used;offset+(modes_used)]),[],3);

% To make the fields "look" continuous for a dielectric we can weight the
% normal field component for the internal field, D=\eps E. \eps\propto k^2
E_int_pm(:,1)=(k_particle/k_medium)^2*E_int_pm(:,1);
E_int_pm2(:,1)=(k_particle/k_medium)^2*E_int_pm2(:,1);

mask=sqrt((X/ac(1)).^2+(Y/ac(1)).^2+(Z/ac(2)).^2)>1;
E_pm=mask(:).*(E_inc+E_sca_pm)+(1-mask(:)).*(E_int_pm);
I_pm=reshape(sum(abs(E_pm).^2,2),size(X));

E_pm2=mask(:).*(E_inc+E_sca_pm2)+(1-mask(:)).*(E_int_pm2);
I_pm2=reshape(sum(abs(E_pm2).^2,2),size(X));


figure(1)
subplot(2,2,2)
imagesc(xv,zv,abs(I_pm))
hold on
plot(majx,majy,'w--');
plot(minx,miny,'w--');
plot(elx,ely,'w--','linewidth',lel);
xlabel('x [\lambda]')
ylabel('z [\lambda]')
grid on
hold off
axis equal square
axis(r_v*[-1,1,-1,1])
caxis(cax)

subplot(2,2,4)
imagesc(xv,zv,abs(I_pm2))
hold on
plot(majx,majy,'w--');
plot(minx,miny,'w--');
plot(elx,ely,'w--','linewidth',lel);
xlabel('x [\lambda]')
ylabel('z [\lambda]')
grid on
hold off
axis equal square
axis(r_v*[-1,1,-1,1])
caxis(cax)

%% radiation patterns
% for the radiation patterns, we also add in the sEBCM and smarties calcs.
% to make things simpler, we will do this in spherical coordinates (doesn't
% make a significant difference in results picking spheroidal coordinates)

ci_in_spherical=find(abs(ab2(1:end/2))|abs(ab2(end/2+1:end)));

[n,m]=combined_index(ci_in_spherical);


T_ebcm=[spheroidal_to_spherical(isProlate,k_medium*c,sT_ebcm(1:end/2,1:end/2),nmax),spheroidal_to_spherical(isProlate,k_medium*c,sT_ebcm(1:end/2,end/2+1:end),nmax);
        spheroidal_to_spherical(isProlate,k_medium*c,sT_ebcm(end/2+1:end,1:end/2),nmax),spheroidal_to_spherical(isProlate,k_medium*c,sT_ebcm(end/2+1:end,end/2+1:end),nmax)];
offset=size(T_sm,1)/2;
T_sma=[T_sm(ci_in_spherical,ci_in_spherical),T_sm(ci_in_spherical,offset+ci_in_spherical);T_sm(offset+ci_in_spherical,ci_in_spherical),T_sm(offset+ci_in_spherical,offset+ci_in_spherical)]


eta=linspace(-1,1,1000)';
xi=inf*ones(size(eta));
phi=zeros(size(eta));
th=acos(eta);


%spheroidal farfield to spherical farfield
[Mff,Nff]=spheroidalvwf_farfield(isProlate,n,m,1e-100,xi,eta,phi,1);

% theta component is -eta component
Mff(:,end/3+1:2*end/3)=-Mff(:,end/3+1:2*end/3);
Nff(:,end/3+1:2*end/3)=-Nff(:,end/3+1:2*end/3);

% vsh are normalised:
Mff=Mff.*1./sqrt(n.*(n+1));
Nff=Nff.*1./sqrt(n.*(n+1));

offset=length(ab2)/2;
abindex=[ci_in_spherical;offset+ci_in_spherical];

% calculate total radiance
I_pm=sum(abs(reshape([Mff.',Nff.']*(ab2(abindex)+2*pq_pm(abindex)),[],3)).^2,2);
I_pm2=sum(abs(reshape([Mff.',Nff.']*(ab2(abindex)+2*pq_pm2(abindex)),[],3)).^2,2);


pq_sm=T_sma*ab2(abindex);
pq_ebcm=T_ebcm*ab2;

I_ebcm=sum(abs(reshape([Mff.',Nff.']*(ab2(abindex)+2*pq_ebcm(abindex)),[],3)).^2,2);
I_sm=sum(abs(reshape([Mff.',Nff.']*(ab2(abindex)+2*pq_sm),[],3)).^2,2);

lntyp={{'k','linewidth',2},...
    {'linestyle',':','linewidth',1},...
    {'linestyle','--','linewidth',1,'color',.5*[1,1,1]},...
    {'linestyle',':','linewidth',1}};

legparams={'sEBCM','sPM','EBCM','PM'};

figure(2)
semilogy(th*180/pi,I_ebcm,lntyp{1}{:})
hold on
semilogy(th*180/pi,I_pm,lntyp{2}{:})
semilogy(th*180/pi,I_sm,lntyp{3}{:})
semilogy(th*180/pi,I_pm2,lntyp{4}{:})
hold off
grid on

xlabel('angle around y-axis');
ylabel('total radiant power [arb]');


I_pm=sum(abs(reshape([Mff.',Nff.']*(2*pq_pm(abindex)),[],3)).^2,2);
I_pm2=sum(abs(reshape([Mff.',Nff.']*(2*pq_pm2(abindex)),[],3)).^2,2);

pq_sm=T_sma*ab2(abindex);
pq_ebcm=T_ebcm*ab2;

I_ebcm=sum(abs(reshape([Mff.',Nff.']*(2*pq_ebcm(abindex)),[],3)).^2,2);
I_sm=sum(abs(reshape([Mff.',Nff.']*(2*pq_sm),[],3)).^2,2);


figure(3)
semilogy(th*180/pi,I_ebcm,lntyp{1}{:})
hold on
semilogy(th*180/pi,I_pm,lntyp{2}{:})
semilogy(th*180/pi,I_sm,lntyp{3}{:})
semilogy(th*180/pi,I_pm2,lntyp{4}{:})
hold off

grid on
xlabel('angle around y-axis');
ylabel('scattered radiant power [arb]');
