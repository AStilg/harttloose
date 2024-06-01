function [M,N,M2,N2,M3,N3] = vswf_vector(n,m,kr,theta,phi,htype)
% vswf.m : Vector spherical wavefunctions: M_k, N_k.
%
% Usage:
% [M,N] = vswf_vector(n,m,kr,theta,phi,type)
% or
% [M1,N1,M2,N2,M3,N3] = vswf_vector(n,m,kr,theta,phi)
% or
% [M1,N1,M2,N2,M3,N3] = vswf_vector(n,kr,theta,phi)
%
% where
% kr, theta, phi are vectors of equal length, or scalar.
% type = 1 -> outgoing solution - h(1)
% type = 2 -> incoming solution - h(2)
% type = 3 -> regular solution - j (ie RgM, RgN)
%
% vector n,m only
%
% M,N are arrays of size length(vector_input,m) x 3
%
% The three components of each vector are [r,theta,phi].
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

% Check input vectors
% These must all be of equal length if non-scalar
% and for good measure, we expand any scalar ones
% to match the others in length

if nargin==3;
    htype=0;
    rtp=kr;
end

if nargin==4
    if ~isscalar(kr)&isscalar(theta)
        htype=theta;
        rtp=kr;
    end
end

if nargin==5
    htype=0;
    rtp=[kr(:),theta(:),phi(:)];
end

if nargin==6
    rtp=[kr(:),theta(:),phi(:)];
end

assert(length(n)==length(m));


Y=zeros(length(m),length(rtp(:,1)));
dtY=Y;
dpY=Y;
hn1=Y;
dhn1=Y;
hn2=Y;
dhn2=Y;
jn=Y;
djn=Y;

for ii=1:max(n)
    findn=find(n==ii);
    if numel(findn)
        mdex=m(findn);
        [Yt,Ythetat,Yphit] = spharm(ii,rtp(:,2),rtp(:,3));

        Y(findn,:)=(Yt(:,ii+(mdex)+1).')*sqrt(ii*(ii+1));
        dtY(findn,:)=(Ythetat(:,ii+(mdex)+1).')/sqrt(ii*(ii+1));
        dpY(findn,:)=(Yphit(:,ii+(mdex)+1).')/sqrt(ii*(ii+1));

        switch htype
            case 1
                [hn1t,dhn1t]=sbesselh1(ii,rtp(:,1));
                hn1(findn,:)=repmat(hn1t.',[length(findn),1]);
                dhn1(findn,:)=repmat(dhn1t.',[length(findn),1]);
            case 2
                [hn1t,dhn1t]=sbesselh2(ii,rtp(:,1));
                hn1(findn,:)=repmat(hn1t.',[length(findn),1]);
                dhn1(findn,:)=repmat(dhn1t.',[length(findn),1]);
            case 3
                [hn1t,dhn1t]=sbesselj(ii,rtp(:,1));
                hn1(findn,:)=repmat(hn1t.',[length(findn),1]);
                dhn1(findn,:)=repmat(dhn1t.',[length(findn),1]);
            otherwise
                [hn1t,dhn1t]=sbesselh1(ii,rtp(:,1));
                hn1(findn,:)=repmat(hn1t.',[length(findn),1]);
                dhn1(findn,:)=repmat(dhn1t.',[length(findn),1]);
                [hn2t,dhn2t]=sbesselh2(ii,rtp(:,1));
                hn2(findn,:)=repmat(hn2t.',[length(findn),1]);
                dhn2(findn,:)=repmat(dhn2t.',[length(findn),1]);
                [jnt,djnt]=sbesselj(ii,rtp(:,1));
                jn(findn,:)=repmat(jnt.',[length(findn),1]);
                djn(findn,:)=repmat(djnt.',[length(findn),1]);
        end

    end
end


%phi elements:
% expimphi=exp(1i*m(:).*rtp(:,3).');
kr(kr==0)=1e-100;
iKR=repmat(1./rtp(:,1).',[length(m),1]);

Y=Y.*iKR;
% dtY=dtY;
% dpY=dpY;


%create the unpacked three vectors for M and N fields:
szj=size(Y,2);

M=zeros(length(n),szj*3);
N=M;

M(:,szj+(1:szj))=hn1.*dpY;
M(:,2*szj+(1:szj))=-hn1.*dtY;

N(:,(1:szj))=hn1.*Y;
N(:,1*szj+(1:szj))=dhn1.*dtY;
N(:,2*szj+(1:szj))=dhn1.*dpY;

M2=0;
N2=0;
M3=0;
N3=0;

if htype==0
    M2=zeros(length(n),szj*3);
    N2=M2;

    M3=N2;
    N3=M3;

    M2(:,szj+(1:szj))=hn2.*dpY;
    M2(:,2*szj+(1:szj))=-hn2.*dtY;

    N2(:,(1:szj))=hn2.*Y;
    N2(:,1*szj+(1:szj))=dhn2.*dtY;
    N2(:,2*szj+(1:szj))=dhn2.*dpY;

    M3(:,szj+(1:szj))=jn.*dpY;
    M3(:,2*szj+(1:szj))=-jn.*dtY;

    N3(:,(1:szj))=jn.*Y;
    N3(:,1*szj+(1:szj))=djn.*dtY;
    N3(:,2*szj+(1:szj))=djn.*dpY;
end

return

