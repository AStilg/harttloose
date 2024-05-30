function F=spheroidal_continued_fractions(isProlate,c,m,lambda,nd,tol)
% spheroidal_continued_fractions - outputs the UNNORMALISED continued
% 								fraction d_{(nd),n}^m(c)/d_{(nd-2),n}^m(c) 
%								correponding to the eigenvalue of spheroidal 
%								mode n,m,c.
%
% USAGE:
%
% F=spheroidal_continued_fractions(isProlate,c,m,lambda,nd,tol);
%
% F 		- d_{(nd),n}^m(c)/d_{(nd-2),n}^m(c) 
% isProlate - true/false
% c         - half interfocal distance.
% m			- azimuthal symmetry.
% lambda    - known eigenvalue
% nd        - target fraction n-mode.
% tol       - set tolerance.
%
% PACKAGE INFO

if nargin<6
    tol=1e-15;
end

sigma=2*isProlate-1;

lambda=lambda+sigma.*c^2; %to put back into offset form.
r=nd-m; %return to flammer form, because the algorithm was confusing me.

%I used mathematica to derive as it took about 2 mins opposed to 10 mins.
% xarg=@(k) (c^4*(-1+2*k-m+nd)*(2*k-m+nd)*(-1+2*k+m+nd)*(2*k+m+nd))/((-3+2*m+2*(2*k-m+nd))*(-1+2*m+2*(2*k-m+nd))^2*(1+2*m+2*(2*k-m+nd)));
xarg=@(k) (sigma.*c^2.*(-1+2*k+2*m+r).*(2*k+2*m+r))./((-1+2*m+2*(2*k+r))*(1+2*m+2*(2*k+r))).*(sigma.*c^2*(-1+2*k+r)*(2*k+r))/((-3+2*m+2*(2*k+r))*(-1+2*m+2*(2*k+r)));
yarg=@(x,k) (2*k+nd)*(1+2*k+nd)+1/2.*sigma.*c^2*(1-(-1+4*m^2)/((3+2*m+2*(2*k-m+nd))*(-1+2*m+2*(2*k-m+nd))))-lambda;

xarg_0=@(k) (sigma.*c^2.*(-1+2*k+r)*(2*k+r))./((-3+2*m+2*(2*k+r))*(-1+2*m+2*(2*k+r)));
if nd==m
    xarg_0=@(k) -(sigma.*c^2./((-3+2*nd)*(-1+2*nd)));
end
if nd==m+1
    xarg_0=@(k) (sigma.*c^2./((-3+2*nd)*(-1+2*nd)));
end


x=1;

[F,err]=continued_fractions_helper(xarg,xarg_0,yarg,x,tol);

if isnan(err)
    F=F(1,:);
    warning('Continued fraction did not converge (is nan)... outputting last result.')
else
    F=F(2,:);
end
% F=F/((c^2*(-1+m+nd)*(m+nd))/((-1+2*m+2*(-m+nd))*(1+2*m+2*(-m+nd)))); % my version of this is included in the recursion... obviously.

end

function [F,err]=continued_fractions_helper(xarg,xarg_0,yarg,x,tol);
x=x(:)';

P=zeros(3,length(x));
Q=zeros(3,length(x));

k=-1;
P(1,:)=1;
Q(1,:)=0;

k=k;
P(2,:)=0*yarg(x,k); %usually this would be the "offset term"... no offset for recursion.
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
