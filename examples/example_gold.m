% example_gold.m
%
% This code was successfully run on a computer with: Intel Core i7 7700K, 32GB ram in MATLAB R2020b.
%
% 

addpath ../ ../provided/ ./helpers/

%distance units in m.
lambda=532e-9;
k_medium=2*pi;
k_particle=2*pi*(.5+2.3i);%(.5+2.3i); 

r=25e-9/lambda; 

AR=linspace(.5,1.5,15);
AR((length(AR)+mod(length(AR),2))/2)=1.00001;
c=0*AR;

% constant volume spheroids
ac=[AR.^(-1/3)*r;AR.^(2/3)*r];


% set nmax
nmax=max(ka2nmax(ac*abs(k_particle)),[],'all')+5;

% generate angular illumination sweep.
th=linspace(0,pi,100)';
ph=0*th;
pol=[1,0].*ones(size(th));

pqCell=cell(length(AR),1);
cdCell=pqCell;
abCell=pqCell;

%generate sT-matrix
warning('off','MATLAB:rankDeficientMatrix');
for ii=1:length(AR)
    [sT_ebcm,sR_ebcm,c(ii)]=stmatrix_spheroid_ebcm(nmax,k_medium,k_particle,ac(:,ii));
    isProlate=AR(ii)>1;
    [a,b]=bsc_plane_spheroidal(isProlate,k_medium*c(ii),nmax,th,ph,pol);

    abCell{ii}=[a;b];
    pqCell{ii}=sT_ebcm*[a;b];
    cdCell{ii}=sR_ebcm*[a;b];

end
warning('on','MATLAB:rankDeficientMatrix');

%% compute scattering cross-section:

% as the spheroidal functions are not properly orthonormal we will have to convert.
C_sca=zeros(length(AR),length(th));
for ii=1:length(AR)

    spq=pqCell{ii};
    sab=abCell{ii};

    nmaxt=floor(sqrt(size(spq,1)/2-1));

    p=zeros(nmaxt*(nmaxt+2),size(spq,2));
    q=p;
    a=p;
    b=q;
    for jj=1:size(spq,2);
        p(:,jj)=spheroidal_to_spherical(AR(ii)>1,k_medium*c(ii),spq(1:end/2,jj),nmaxt);
        q(:,jj)=spheroidal_to_spherical(AR(ii)>1,k_medium*c(ii),spq(end/2+1:end,jj),nmaxt);
        a(:,jj)=spheroidal_to_spherical(AR(ii)>1,k_medium*c(ii),sab(1:end/2,jj),nmaxt);
        b(:,jj)=spheroidal_to_spherical(AR(ii)>1,k_medium*c(ii),sab(end/2+1:end,jj),nmaxt);
    end
    C_sca(ii,:)=sum(abs(p).^2+abs(q).^2)./sum(abs(a).^2+abs(b).^2);
end
%%

%set the plots and beam angle
plot_field_index=[1,8,15];
beam_index=10;
th(beam_index)*180/pi;

[ARG,TH]=meshgrid(AR,th);
figure(100)
set(100,'position',[285   283   809   276])
subplot(1,3,1)
surf(TH*180/pi,ARG,C_sca.','EdgeColor','none')
hold on
plot3(TH(:,plot_field_index)*180/pi,ARG(:,plot_field_index),C_sca(plot_field_index,:).','k','LineWidth',2)
plot3(TH(beam_index,:).'*180/pi,ARG(beam_index,:),C_sca(:,beam_index),'r','LineWidth',2)
hold off
axis square
grid on
xlim([0,180])
xlabel('incident planewave angle [deg]')
ylabel('aspect ratio');
zlabel('C_{sca}');
xticks([0:30:180]);
subplot(1,3,2)
surf(TH*180/pi,ARG,C_sca.',TH*180/pi,'edgecolor','none')
hold on
plot3(TH(:,plot_field_index)*180/pi,ARG(:,plot_field_index),C_sca(plot_field_index,:).','k','LineWidth',2)
% plot3(TH(beam_index,:).'*180/pi,ARG(beam_index,:),C_sca(:,beam_index),'r','LineWidth',2)
hold off
axis square
grid on
xlim([0,180])
xlabel('incident planewave angle [deg]')
ylabel('aspect ratio');
zlabel('C_{sca}');
xticks([0:30:180]);
view(90,0)
subplot(1,3,3)
surf(TH*180/pi,ARG,C_sca.',ARG,'edgecolor','none')
hold on
plot3(TH(:,plot_field_index)*180/pi,ARG(:,plot_field_index),C_sca(plot_field_index,:).','k','LineWidth',2)
plot3(TH(beam_index,:).'*180/pi,ARG(beam_index,:),C_sca(:,beam_index),'r','LineWidth',2)
hold off
axis square
grid on
xlim([0,180])
xlabel('incident planewave angle [deg]')
ylabel('aspect ratio');
zlabel('C_{sca}');
xticks([0:30:180]);
view(0,0)
% exportgraphics(gcf,'Csca.png','resolution',300)

%% check field continuity at extremum.

%get range for xi.

%we want a square region around the ellipsoid. first find biggest
%dimension:

rmax=1.05*max(ac,[],'all');

thcirc=linspace(0,2*pi);

[n,m]=combined_index([0:nmax*(nmax+2)]');


figure(10)
set(10,'position',[285   283   809   276])
for jj=1:length(plot_field_index)
    ii=plot_field_index(jj);

    spq=pqCell{ii};
    sab=abCell{ii};
    scd=cdCell{ii};

    isProlate=AR(ii)>1;
    
    [extMesh,intMesh]=nice_spheroidal_mesh(ac(:,ii),21,20,2*1.05*max(ac(:)),2*1.05*max(ac(:)));

    [Msca,Nsca,~,~,Minc,Ninc]=spheroidalvwf(isProlate,n,m,k_medium*c(ii),extMesh.XI(:),extMesh.ETA(:),extMesh.PHI(:));
    [Mint,Nint]=spheroidalvwf(isProlate,n,m,conj(k_particle)*c(ii),intMesh.XI(:),intMesh.ETA(:),intMesh.PHI(:),3); %why conjugate?

    H_ext=reshape(1*[Nsca.',Msca.']*spq(:,beam_index)+1*[Ninc.',Minc.']*sab(:,beam_index),[],3);
    H_int=(k_particle/k_medium)^1*reshape(1*[Nint.',Mint.']*scd(:,beam_index),[],3);

    I_H_ext=(sum(abs(H_ext).^2,2));
    I_H_int=(sum(abs(H_int).^2,2));
    
    subplot(1,length(plot_field_index),jj)
    surf(extMesh.X,extMesh.Z,extMesh.Y,(reshape(I_H_ext,size(extMesh.XI))),'edgecolor','interp','facecolor','interp')
    hold on
    surf(intMesh.X,intMesh.Z,intMesh.Y,(reshape(I_H_int,size(extMesh.XI))),'edgecolor','interp','facecolor','interp')
    plot(ac(1,ii)*cos(thcirc),ac(2,ii)*sin(thcirc),'w--')
    hold off
    view(2)
    axis equal
    xlim([-1.05*max(ac(:)),1.05*max(ac(:))])
    ylim([-1.05*max(ac(:)),1.05*max(ac(:))])
    grid on
    xlabel('x-position [\lambda]')
    ylabel('z-position [\lambda]')
end
% exportgraphics(gcf,'goldnf.png','Resolution',300)
