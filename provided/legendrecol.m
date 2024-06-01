function pnm=legendrecol(n,theta,m,xflag)
% legendrecol.m : Gives the spherical coordinate recursion in n for a given
%                 m, theta. m=0 unless specified.
%
% Usage:
% pnm = legendrecol(n,theta,m); %spherical harmonic legendre.
% pnm = legendrecol(n,theta,m,1); %theta->x, |x|>0 legendre.
%
% Inspiration from [Holmes and Featherstone, 2002]. Note: row recursions
% limit to [0:m], col recursions all m=m.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

sgnm=1;
if (m<0)
    sgnm=(-1).^m;
end
m=abs(m);

if n==0;
    pnm=1/sqrt(2*pi)/sqrt(2);
    return %don't waste my time.
end

if nargin<4
    xflag=0;
    if nargin<3
        m=0;
    end
end

theta=theta(:).';

ct=cos(theta);
st=sin(theta);
if xflag
    ct=theta;
    [N,D]=rat(real(theta),1e-15);
    xp=(D+N)./D+1i*imag(theta);
    xm=(D-N)./D-1i*imag(theta);
    st=sqrt(xp.*xm);
    % alternate way is to swap theta and 1.
end

%column starts along the m=n diagonal:
Wmm=sqrt((2*m+1)/(4*pi)*prod(1-1/2./[1:m]))*ones(size(theta)); %first entry! %*sin^m omitted as it's a "common" factor.
Wmp1m=sqrt(2*m+3).*ct.*Wmm; %relationship between mm and mp1m terms.

lnm=length([m:n]); % we're going to n from m.

pnm=zeros(lnm,length(theta));

pnm(1,:)=Wmm;

%recursion occurs down the column:
if lnm~=1
    pnm(2,:)=Wmp1m;
    for ii=2:n-m %compute to n including m, ii is current n-m.
        jj=ii+m; %actual n index with m.
        a=sqrt((2*jj-1).*(2.*jj+1)./(jj-m)/(jj+m));
        b=sqrt((2*jj+1)*(jj-m-1)*(jj+m-1)./(jj-m)/(jj+m)/(2*jj-3));
        pnm(ii+1,:)=a*ct.*pnm(ii,:)-b*pnm(ii-1,:);
    end
end

pnm=sgnm.*pnm.*st.^m; %re-introduce sin(theta).^m... for reasons. I.e. no idea why it has to be done this way (it should work starting from oscillating anyway?).
% pnm=pnm(end,:);