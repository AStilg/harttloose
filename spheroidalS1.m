function [Snm,dSnm,dphiSnm,eigenvalues,Norm]=spheroidalS1(isProlate,n,m,c,eta,allvalues);
% spheroidalS1 - calculates the normalised angular spheroidal
%				wavefunctions, \bar S_{nm}(c,\eta). Can output all
%				wavefunctions of same m and partiy.
%
% usage:
% [Snm,dSnm,dphiSnm,eigenvalues]=spheroidalS1(isProlate,n,m,c,eta)
%
% or 
% 
% [Snm,dSnm,dphiSnm,eigenvalues]=spheroidalS1(isProlate,n,m,c,eta,allvalues)
%
% Default behaviour is that only the requested value is given.
%
% NOTE: The derivative is,
%
% dSnm = \sqrt{1-\eta^2}\frac{d}{d\eta} \bar S_{nm}(c,\eta)
%  and
% dphiSnm = m*Snm/\sqrt{1-x.^2}
% m>=0 only
%
% NOTE: if allvalues, the last few rows are guaranteed to be wrong.
%
% PACKAGE INFO

debug=false;

if nargin<6
    allvalues=false;
end

if lower(allvalues)=='all'
    allvalues=true;
end

if n<abs(m)
    error('n must be greater or equal to |m|.')
end

if m<0
    error('m must be greater or equal to 0. dael with signs in derived functions.')
end

% derived parameters:
eta=eta(:).';
isodd=mod(n+abs(m),2);
N=ceil(n+abs(c)+15);

% compute superposition weights:
if c~=0
    [U,r,eigenvalues]=spheroidal_u_coefficients(isProlate,isodd,abs(m),c,N);
else
    U=eye(N+1);
    r=2*[0:N]'+isodd+abs(m);
    eigenvalues=r.*(r+1);
end

if debug
    figure
    imagesc(r',r,angle(U))
end

if ~allvalues
    U=(U((r==n),:));
    eigenvalues=eigenvalues((r==n));
else
    
%     n=r';
%     if (m==0)
%         U(:,1)=0;
%     end
end

%create storage for associated legendre functions:
temp=legendrecol(r(end)+1+2,eta,abs(m),1);
Y=temp(r-abs(m)+1,:);

%compute Snm
Snm=(-1).^(m).*(U*Y); %good.

if nargout<2
    return;
end

%generate the reference values:
tempp=legendrecol(r(end)+1+2,eta,abs(m)+1,1); % will output one lessor n0=m+1
tempm=legendrecol(r(end)+1+2,eta,abs(m)-1,1); % will output one more n0=m-1 m=0 is special case where both lessor.

%calculate the n-vectors used.
nv=[abs(m):r(end)+1]';
nv=nv(r-abs(m)+1);

% compute -1i*phi derivative if needed?
dpY=(1/2 * sqrt((2*nv+1)./(2*nv+3)) .* ...
    ( sqrt((nv+abs(m)+1).*(nv+abs(m)+2)) .* tempp(r-abs(m)+1,:) ...
    + sqrt((nv-abs(m)+1).*(nv-abs(m)+2)) .* tempm(r-abs(m)+1+2*(m~=0),:) ));

% compute theta derivative. Padding to so I can re-use index.
tempp=[0*tempp(1,:);tempp]; % this will always be padded.
if m==0
    tempm=[0*tempm(1,:);tempm]; % this is padded if m=0.
end
dtY =-(sqrt((nv-abs(m)+1).*(nv+abs(m)))/2 .* tempm(r-abs(m)+2-(m==0),:) ...
    - sqrt((nv-abs(m)).*(nv+abs(m)+1))/2 .* tempp(r-abs(m)+1,:));

%compute dSnm
dSnm=(-1).^(m).*(U*dtY); %good.

%compute dphiSnm
dphiSnm=(-1).^(m).*(U*dpY); %good.

if nargout==5
    temp=abs(U).^2.*r'.*(r'+1);
    if m==0
        temp(:,r==0)=abs(U(:,r==0)).^2.*sqrt(4*pi);
    end
    Norm=sum(temp,2);
end


