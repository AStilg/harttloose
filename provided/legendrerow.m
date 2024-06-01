function pnm=legendrerow(n,theta,m,xflag)
% legendrerow.m : Gives the spherical coordinate recursion in m for a given
%                 n, theta.
%
% Usage:
% pnm = legendrerow(n,theta)
%
% This provides approximately no benefit over the MATLAB implimentation. It
% *may* provide a benefit in Octave. Inspiration from [Holmes and Featherstone, 2002]
% and [Jekeli et al., 2007].
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if n==0;
    pnm=1/sqrt(2*pi)/sqrt(2);
    return
end

if nargin<4
    xflag=0;
    if nargin<3
        m=n;
    end
end

theta=theta(:).'; %bug fix for complex case.

ct=cos(theta);
st=sin(theta);
if xflag
    ct=theta;
    
    [N,D]=rat(real(theta),1e-15);
    xp=(D+N)./D+1i*imag(theta);
    xm=(D-N)./D-1i*imag(theta);
    st=sqrt(xp.*xm);

end

Wnn=sqrt((2*n+1)/(4*pi)*prod(1-1/2./[1:n]))*ones(size(theta)); %first entry!
Wnnm1=sqrt(2*n)*ct.*Wnn; %second entry!
lnm=length([0:n]);

pnm=zeros(lnm,length(theta));
pnm(end,:)=Wnn;
pnm(end-1,:)=Wnnm1;

if lnm==2;
    pnm=[Wnnm1;Wnn];
else
    jj=lnm-2;
    for ii=n-2:-1:0
        a=sqrt(4*(ii+1)^2/(n-ii)/(n+ii+1));
        b=sqrt((n-ii-1)*(n+ii+2)/(n-ii)/(n+ii+1));
        
        pnm(jj,:)=a*ct.*pnm(jj+1,:)-b*st.^2.*pnm(jj+2,:);   %row recursion!
        jj=jj-1;
    end
end

[ST,M]=meshgrid(st,(0:m));

pnm=pnm(1:m+1,:).*ST.^(M);


