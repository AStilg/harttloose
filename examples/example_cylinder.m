% example_cylinder.m
%
% This code was successfully run on a computer with: Intel Core i7 7700K, 32GB ram in MATLAB R2020b.
%
% 


addpath ../ ../provided/ ./helpers/

%all distance units are \lambda.
k_medium=2*pi;
k_particle=2*pi*1.5;

% set up width of cylinder
r0=1;
z0=r0*logspace(-1,1,5);
z0(z0==r0)=z0(z0==r0)*1.00001; %so that the cylinder is not precisely 1/1 AR.
wh=[r0*ones(size(z0));z0];

%calculate an estimate for c .
[isP,cref]=aspect_ratio_to_conk(wh/2); %is done in terms of radius not diameters.

T1=cell(size(wh,2),1);
T2=T1;
% T3=T2;
sT1=cell(size(wh,2),1);
sT2=sT1;

c1=zeros(size(wh,2),1);
c2=zeros(size(wh,2),1);

warning('off','MATLAB:rankDeficientMatrix')
nmaxoffset=[0,0,0,1,3];
for ii=1:size(wh,2)
    ii
    nmax1=ka2nmax(pi*max(wh(:,ii)));
    snmax=nmax1-nmaxoffset(ii);

    isProlate=wh(2,ii)>wh(1,ii);

    [T1{ii}]=tmatrix_cylinder_ebcm_romberg(nmax1,k_medium,k_particle,wh(:,ii));
    [T2{ii}]=tmatrix_cylinder_pm(nmax1,k_medium,k_particle,wh(:,ii));
    [sT1{ii},~,c1(ii)]=stmatrix_cylinder_ebcm(snmax,k_medium,k_particle,wh(:,ii));
    [sT2{ii},~,c2(ii)]=stmatrix_cylinder_pm(snmax,k_medium,k_particle,wh(:,ii));
end
warning('on','MATLAB:rankDeficientMatrix')

%% Set up beam for spherical we are going to use spheroidal functions in the far-field limit. This will be fine.

angl=pi/4;
pol=[1,0];
l=0;

nmax1=max(ka2nmax(pi*max(wh)));
[a,b]=bsc_bessel_spheroidal(isProlate,0,nmax1,angl,l,pol);
%normalise to vswf:
a=spheroidal_to_spherical(isProlate,0,a,nmax1);
b=spheroidal_to_spherical(isProlate,0,b,nmax1);

eta=linspace(-1,1,1000)';
xi=inf*ones(size(eta));
phi=zeros(size(eta));
th=acos(eta);

ci_in_spherical=find(abs(a)|abs(b));
[n,m]=combined_index(ci_in_spherical);

%spheroidal farfield to spherical farfield
[Mff,Nff]=spheroidalvwf_farfield(isProlate,n,m,1e-100,xi,eta,phi,1);

% theta component is -eta component
Mff(:,end/3+1:2*end/3)=-Mff(:,end/3+1:2*end/3);
Nff(:,end/3+1:2*end/3)=-Nff(:,end/3+1:2*end/3);

%normalise:
Mff=Mff.*1./sqrt(n.*(n+1));
Nff=Nff.*1./sqrt(n.*(n+1));


lntyp={{'k','linewidth',2},...
    {'linestyle',':','linewidth',1},...
    {'linestyle','--','linewidth',1,'color',.5*[1,1,1]},...
    {'linestyle',':','linewidth',1}};

legparams={'sEBCM','sPM','EBCM','PM'};


figure(2020)
clf
for ii=1:size(wh,2)
    %% spherical part.
    nmax1=floor(sqrt(size(T1{ii})/2));

    [n,m]=meshgrid([min([abs(l),1]):nmax1],l+[-1,1]);
    n=n(:);
    m=m(:);

    remel=find(n<abs(m));

    n(remel)=[];
    m(remel)=[];

    indexing=combined_index(n(:),m(:));

    pq_ebcm=[T1{ii}(indexing,indexing),T1{ii}(indexing,indexing+end/2);T1{ii}(indexing+end/2,indexing),T1{ii}(indexing+end/2,indexing+end/2)]*[a(indexing);b(indexing)];
    pq_pm=[T2{ii}(indexing,indexing),T2{ii}(indexing,indexing+end/2);T2{ii}(indexing+end/2,indexing),T2{ii}(indexing+end/2,indexing+end/2)]*[a(indexing);b(indexing)];

    I_sca_ebcm=sum(abs(reshape([Mff(1:length(pq_pm)/2,:).',Nff(1:length(pq_pm)/2,:).']*pq_ebcm,[],3)).^2,2);
    I_sca_pm=sum(abs(reshape([Mff(1:length(pq_pm)/2,:).',Nff(1:length(pq_pm)/2,:).']*pq_pm,[],3)).^2,2);
    %% spheroidal part
    snmax1=floor(sqrt(size(sT1{ii})/2-1));

    [sn,sm]=meshgrid([0:snmax1],l+[-1,1]);
    sn=sn(:);
    sm=sm(:);

    remel=find(sn<abs(sm));

    sn(remel)=[];
    sm(remel)=[];
    indexing=combined_index(sn(:),sm(:))+1;

    isProlate=wh(2,ii)>wh(1,ii);
    
    [sMff,sNff]=spheroidalvwf_farfield(isProlate,sn,sm,k_medium*c1(ii),xi,eta,phi,1);

    [sa,sb]=bsc_bessel_spheroidal(isProlate,k_medium*c1(ii),snmax1,angl,l,pol);

    spq_ebcm=[sT1{ii}(indexing,indexing),sT1{ii}(indexing,indexing+end/2);sT1{ii}(indexing+end/2,indexing),sT1{ii}(indexing+end/2,indexing+end/2)]*[sa(indexing);sb(indexing)];
    spq_pm=[sT2{ii}(indexing,indexing),sT2{ii}(indexing,indexing+end/2);sT2{ii}(indexing+end/2,indexing),sT2{ii}(indexing+end/2,indexing+end/2)]*[sa(indexing);sb(indexing)];

    sI_sca_ebcm=sum(abs(reshape([sMff(1:length(spq_pm)/2,:).',sNff(1:length(spq_pm)/2,:).']*spq_ebcm,[],3)).^2,2);
    sI_sca_pm=sum(abs(reshape([sMff(1:length(spq_pm)/2,:).',sNff(1:length(spq_pm)/2,:).']*spq_pm,[],3)).^2,2);

    subplot(2,3,ii)
    semilogy(th*180/pi,sI_sca_ebcm,lntyp{1}{:})
    hold on
    semilogy(th*180/pi,sI_sca_pm,lntyp{2}{:})
    semilogy(th*180/pi,I_sca_ebcm,lntyp{3}{:})
    semilogy(th*180/pi,I_sca_pm,lntyp{4}{:})
    yl=ylim;
    grid on
    val=-4-1*(ii==5);
    ylim([10^val,1]*10^(ceil(log10(yl(2)))))
    xlim([0,180])
    yticks([10.^([ceil(val/2)*2:2:0])*10^ceil(log10(yl(2)))])
    plot(180-45*[1,1],ylim,'color',[.4,.4,.4])
    hold off
    % ylim([1e-5,10000])
    xlabel('angle around y-axis');
    ylabel('scattered radiant power [arb]');
end

subplot(2,3,6)
plot(0,0,lntyp{1}{:})
hold on
plot(0,0,lntyp{2}{:})
plot(0,0,lntyp{3}{:})
plot(0,0,lntyp{4}{:})
hold off
h=legend(legparams{:});
axis off

% exportgraphics(gcf,'cylinderscattering.pdf')