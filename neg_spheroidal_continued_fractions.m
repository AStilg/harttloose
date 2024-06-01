function F=neg_spheroidal_continued_fractions(isProlate,c,m,lambda,nd,tol)
% neg_spheroidal_continued_fractions - outputs the UNNORMALISED continued
% 								fraction d_{(-nd),n}^m(c)/d_{(-nd+2),n}^m(c) 
%								correponding to the eigenvalue of spheroidal 
%								mode n,m,c.
% note d_{-nd,n}^{m}/d_{-nd+2,n}^{m} -> d_{-r,n}^{m}/d_{-r+2,n}^{m} in Flammer.
%
% USAGE:
%
% F=neg_spheroidal_continued_fractions(isProlate,c,m,lambda,nd,tol);
%
% F 		- d_{(-nd),n}^m(c)/d_{(-nd+2),n}^m(c)
% isProlate - true/false
% c         - half interfocal distance.
% m			- azimuthal symmetry.
% lambda    - known eigenvalue
% nd        - target fraction n-mode.
% tol       - set tolerance.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if nargin<6
    tol=1e-15;
end

sigma_c=2*isProlate-1;

lambda=lambda+sigma_c.*c^2; %to put back into offset form.

r=nd-abs(m); %algorithm was confusing me so I re-defined to match flammer.

sigma=mod(r+1,2);
sgn=(1-2*sigma);

%mathematica to the rescue.
xarg=@(k) (sigma_c.*c^2*(1-2*k+2*m-r)*(2-2*k+2*m-r))/((-1+2*m+2*(2-2*k-r))*(1+2*m+2*(2-2*k-r)))*(sigma_c.*c^2*(1-2*k-r)*(2-2*k-r))/((1+2*m+2*(-2*k-r))*(3+2*m+2*(-2*k-r))); %A*C
yarg=@(x,k) -lambda+(sigma_c.*c^2*(-1-2*m^2+2*(-2*k+m-r)*(1-2*k+m-r)))/((-1+2*m+2*(-2*k-r))*(3+2*m+2*(-2*k-r)))+(-2*k+m-r)*(1-2*k+m-r); %B

xarg_0=@(k) (sigma_c.*c^2*(1-2*k+2*m-r)*(2-2*k+2*m-r))/((-1+2*m+2*(2-2*k-r))*(1+2*m+2*(2-2*k-r))); %A

if r==2*m+1+sigma
%     disp('special case')
%     xarg_0=@(k) sgn*c^2/((r-3*sigma+2*k+2).*(r-2+2*k-3*sigma+2));
    xarg_0=@(k) sigma_c.*sgn*c^2/((r+1*sigma+2*k-2).*(r-4+2*k+1*sigma)); %A otherwise for problematic case.
end
x=1;

[F,err]=continued_fractions_helper_neg(xarg,xarg_0,yarg,x,tol);

if isnan(err)
    F=F(1,:);
    warning('Continued fraction did not converge (is nan)... outputting last result.')
else
    F=F(2,:);
end

end

function [F,err]=continued_fractions_helper_neg(xarg,xarg_0,yarg,x,tol);
x=x(:)';

P=zeros(3,length(x));
Q=zeros(3,length(x));

k=-1;
P(1,:)=1;
Q(1,:)=0;

k=k;
P(2,:)=0; %usually this would be the "offset term"... no offset for recursion.
Q(2,:)=1;

k=k+1;
P(3,:)=yarg(x,k).*P(2,:)-xarg_0(k)*P(1,:);
Q(3,:)=yarg(x,k).*Q(2,:)-xarg_0(k)*Q(1,:);

F=P(2:3,:)./Q(2:3,:);
go=any(abs(F(2,:)-F(1,:))>abs(F(2,:))*tol);

while and(go,~any(isnan(F(:))))

    P(1:2,:)=P(2:3,:);
    Q(1:2,:)=Q(2:3,:);
    k=k+1;
    P(3,:)=(yarg(x,k).*P(2,:)-xarg(k)*P(1,:)); %negative ratio
    Q(3,:)=(yarg(x,k).*Q(2,:)-xarg(k)*Q(1,:));
    F=P(2:3,:)./Q(2:3,:);
    
    go=any(abs(F(2,:)-F(1,:))>abs(F(2,:))*tol);
end

%error
err=abs(F(2,:)-F(1,:))/abs(F(2,:));

end
