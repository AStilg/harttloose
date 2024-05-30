function [Y,Ytheta,Yphi] = spharm(n,m,theta,phi,varargin)
% spharm.m : scalar spherical harmonics and
%            angular partial derivatives for given n,m (can take vector m).
%
% Usage:
% [Y,dY/dtheta,1/sin(theta)*dY/dphi] = spharm(n,m,theta,phi)
% or
% [Y,dY/dtheta,1/sin(theta)*dY/dphi] = spharm(n,theta,phi)
% or
% [Y,dY/dtheta,1/sin(theta)*dY/dphi] = spharm(n,theta);
%
% Scalar n for the moment.
%
% If scalar m is used Y is a vector of length(theta,phi) and is
% completely compatible with previous versions of the toolbox. If vector m
% is present the output will be a matrix with rows of length(theta,phi) for
% m columns.
%
% "Out of range" n and m result in return of Y = 0
%
% PACKAGE INFO

% normalise=0;
if length(n)>1
    error('n must be a scalar.')
end

if nargin==2
    phi=0;
    theta=m;
    m=[-n:n];
end

if nargin==3
    phi=theta;
    theta=m;
    m=[-n:n];
end

%this is a cop out meant for future versions.
if nargout>1
    mi=m;
    m=[-n:n];
end

m=m(abs(m)<=n);

[theta,phi] = matchsize(theta,phi);

% if abs(m) > n | n < 0
%    Y = zeros(input_length,1);
%    Ytheta = zeros(input_length,1);
%    Yphi = zeros(input_length,1);
%    return
% end
if n>0
    pnm = legendrerow(n,theta);
else
    pnm = ones(size(theta))/sqrt(4*pi);
end

pnm = pnm(abs(m)+1,:); %pick the m's we potentially have. 

[phiM,mv]=meshgrid(phi,m);

expphi = exp(1i*mv.*phiM);

pnm = [(-1).^mv(m<0,:).*pnm(m<0,:);pnm(m>=0,:)].*(-1).^(mv);%.*(-1).^(mv)

Y = pnm .* expphi;

%N = sqrt((2*n+1)/(8*pi));


% Do we want to calculate the derivatives?
if nargout <= 1
    Y=Y.';
%     if normalise
%         Y=Y/sqrt(n*(n+1));
%     end
    % Doesn't look like it
    return
end

% We use recursion relations to find the derivatives, choosing
% ones that don't involve division by sin or cos, so no need to
% special cases to avoid division by zero

% exp(i*phi),exp(-i*phi) are useful for both partial derivatives
expplus = exp(1i*phiM);
expminus = exp(-1i*phiM);

% theta derivative
% d/dtheta Y(n,m) = 1/2 exp(-i phi) sqrt((n-m)(n+m+1)) Y(n,m+1)
%                 - 1/2 exp(i phi) sqrt((n-m+1)(n+m)) Y(n,m-1)

ymplus=[Y(2:end,:);zeros(1,length(theta))];
ymminus=[zeros(1,length(theta));Y(1:end-1,:)];

Ytheta = -(sqrt((n-mv+1).*(n+mv))/2 .* expplus .* ymminus ...
    - sqrt((n-mv).*(n+mv+1))/2 .* expminus .* ymplus);

% phi derivative - actually 1/sin(theta) * d/dphi Y(n,m)
% Note that this is just i*m/sin(theta) * Y(n,m), but we use a
% recursion relation to avoid divide-by-zero trauma.
% 1/sin(theta) d/dphi Y(n,m) =
% i/2 * [ exp(-i phi) sqrt((2n+1)(n+m+1)(n+m+2)/(2n+3)) Y(n+1,m+1)
%     + exp(i phi) sqrt((2n+1)(n-m+1)(n-m+2)/(2n+3)) Y(n+1,m-1) ]

Y2 = spharm(n+1,theta,phi).';

ymplus=Y2(3:end,:);
ymminus=Y2(1:end-2,:);

Yphi = -(1i/2 * sqrt((2*n+1)/(2*n+3)) * ...
    ( sqrt((n+mv+1).*(n+mv+2)) .* expminus .* ymplus ...
    + sqrt((n-mv+1).*(n-mv+2)) .* expplus .* ymminus ));

Y=(Y(n+mi+1,:)).';
Yphi=(Yphi(n+mi+1,:)).';
Ytheta=(Ytheta(n+mi+1,:)).';

% if normalise
%     Y=Y/sqrt(n*(n+1));
%     Yphi=Yphi/sqrt(n*(n+1));
%     Ytheta=Ytheta/sqrt(n*(n+1));
% end

return
