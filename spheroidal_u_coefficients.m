function [u,n,eigenvalues,W,uold] = spheroidal_u_coefficients(isProlate,isOdd,m,c,N);
% spheroidal_u_coefficients - computed the orthonormal superposition of
% 							spherical harmonics to generate a normalised
%						    S_nm(c,eta).
%
% usage:
%
% [U,n,eigenvalues] = spheroidal_u_coefficients(isProlate,isOdd,m,c,N);
%
% The ellipticity factor is "c", "m" is the azimuthal mode, "isodd" is the
% parity of mod(n+m,2), "N" is the number of modes needed---if blank, it
% uses N=ceil(abs(c)+50).
%
% Usage of outputs:
%
% \bar{S}_n(i)|m|(c,eta)=\sum_i^N u(i,j)*Y_n(j)|m|(eta,0).
%
% where \bar{S}_n(i)m is the normalised angular spheroidal harmonic, Y_n(j)m is
% the normalised spherical harmonic. Note: \bar{S}_n(j)m will obtain a
% factor of 1/sqrt(2*pi) due to the normalisation of Y_im.
%
% Parity switches between either odd or even series in n. Thus, the columns
% of u in terms of mode n are: n(i)\in[m+isodd:2:m+2*N+isodd].
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if nargin<4
    N=ceil(abs(c)+35);
end

n=[0:N]';
r=2*n+isOdd;
m=abs(m);

sigma=2*isProlate-1;

%create n indexes
n=r+abs(m);

if abs(c)~=0
    %compute coupling coefficients for recursion:
%     A=(2*m+r+2).*(2*m+r+1)./((2*m+2*r+3).*(2*m+2*r+5));
%     B=(m+r).*(m+r+1).*sigma./c.^2+(2*(m+r).*(m+r+1)-2*m.^2-1)./((2*m+2*r-1).*(2*m+2*r+3))-1; %we've included an offset. It may help?
%     C=r.*(r-1)./((2*m+2*r-3).*(2*m+2*r-1));

    %symmetrize:
%     AC=((1+r).*(2+r).*(1+2*m+r).*(2+2*m+r))./((1+2*m+2.*r).*(3+2*m+2.*r).^2.*(5+2*m+2.*r));
%     dm=sqrt(AC(1:end-1));
    
    %construct eigenproblem
%     EigenProblem=diag(B)+diag(dm,1)+diag(dm,-1);

    %normalised values
%     A=(2*m+r+2).*(2*m+r+1)./((2*m+2*r+3).*(2*m+2*r+5)).*sqrt(((-2+m-n).*(-1+m-n).*(5+2.*n))./((1+m+n).*(2+m+n).*(1+2.*n)));
    A=sqrt((1+r).*(2+r).*(2*m+r+1).*(2*m+r+2)./(1+2.*(m+r))./(2*m+2*r+5))./(2*m+2*r+3); %let's simplify!
    
    B=((m+r).*(m+r+1)./c.^2.*sigma+(2*(m+r).*(m+r+1)-2*m.^2-1)./((2*m+2*r-1).*(2*m+2*r+3))-1);
%     C=r.*(r-1)./((2*m+2*r-3).*(2*m+2*r-1)).*sqrt(((-1+m+n).*(m+n).*(-3+2.*n))./((m-n).*(1+m-n).*(1+2.*n))); 
    
    EigenProblem=diag(B)+diag(A(1:end-1),1)+diag(A(1:end-1),-1);
    
    %solve eigensystem
    [u,eigenvalues]=eig(EigenProblem,'nobalance');
    
    ev0=eigenvalues;
    
    eigenvalues=diag(eigenvalues).*sigma.*c.^2;
    
    %let's sort eigenvalues
    [~,sorted_index]=sort(real(eigenvalues));
    eigenvalues=eigenvalues(sorted_index);
    ev0=ev0(sorted_index,:);
    u=u(:,sorted_index);
    
    uold=u.';
    
    % now refine the elements, must be done in MATLAB (but not octave) because
    % it uses an approximate method to find eigenvectors which produces LARGE
    % errors in integrals, etc. tested in R2020b.
    u_rat_down=zeros(1,size(u,2)-1);
    u_rat_up=zeros(1,size(u,2)-1);
    for ii=1:length(u_rat_down)
        u_rat_down(ii)=spheroidal_continued_fractions(isProlate,c,m,eigenvalues(ii),n(end))./sqrt(((m-n(end)).*(1+m-n(end)).*(1+2.*n(end)))./((-1+m+n(end)).*(m+n(end)).*(-3+2.*n(end))));
        %     pause
        u_rat_up(ii)=sqrt(((-2+m-n(1)).*(-1+m-n(1)).*(5+2.*n(1)))./((1+m+n(1)).*(2+m+n(1)).*(1+2.*n(1)))).*neg_spheroidal_continued_fractions(isProlate,c,m,eigenvalues(ii+1),-n(1)+2*m);
        %     pause
    end
    
    u(1,2:end)=u_rat_up;
    for ii=2:size(u,2)-1
        u(ii,ii+1:end)=up_recursion_ratio_norm(sigma,n(ii),m,c,eigenvalues(ii+1:end).',u(ii-1,ii+1:end));
    end
    
    u(end,1:end-1)=u_rat_down;
    for ii=size(u,2)-1:-1:2
        u(ii,1:ii-1)=down_recursion_ratio_norm(sigma,n(ii),m,c,eigenvalues(1:ii-1).',u(ii+1,1:ii-1));
    end
    
    for jj=1:size(u,2)-1
        u(jj:end,jj)=cumprod(u(jj:end,jj));
        u(jj+1:-1:1,jj+1)=cumprod(u(jj+1:-1:1,jj+1));
    end
    
    % % check eigenvectors, one-by-one:
    % for ii=1:size(u,2)
    %    mean([u(:,ii)./(EigenProblem*u(:,ii))*eigenvalues(ii)])
    % end    
    
    u=u.'; %transpose so that matrix multiplication may be applied in calcs.
    
    %hack the coefficients using the signed real parts of the coefficient,
    %this ensures that there is consistency between different problem
    %sizes. This is not without consequence. It means that problems may
    %spontaneously change sign. A better fix may be to use extremum values
    %of the angular spheroidal wavefunctions. A good choice for this will
    %be eta=0 for even functions and gradients for odd.    
    
	Pmm=(-1)^(m).*sqrt((2*m+1)./(4*pi)).*prod([1:2:2*m-1])./prod(sqrt([1:2*m])); %calculate ass legendre function P_m^m(0).

    v=[Pmm;-sqrt(((-1+m-(n-isOdd)).*(1+m+(n-isOdd)).*(5+2.*(n-isOdd)))./((-2+m-(n-isOdd)).*(2+m+(n-isOdd)).*(1+2.*(n-isOdd))))];
	
	cv=cumprod(v);
    
    if isOdd
        v=-sqrt(m+n).*sqrt((-m+n).*(1+2.*n)./((-1+2.*n))).*cv(1:end-1);
    else
        v=cv(1:end-1);
    end
    
    %generate x=0 to correct sign;
    x=u*v;
    
    signs=sign(real(x).*real(v));
    signs(signs==0)=1;
    u=u.*signs; 
    
    %override to unity middle:
%     u=u.*sign(real(diag(u)));
    
else
    u=eye(length(n));
    eigenvalues=n.*(n+1);
end

% output the matrix that can yield d=W.*u.
if nargout>3
    tt=spherical_harmonic_normalisation(n(:),m);
    W=tt(1:end)'./tt(1:end);
end

end


function u1u0=down_recursion_ratio_norm(sigma,nd,m,c,ev,u2u1);
% note that this is u_nd/u_nd-2. input is u_nd+2/u_nd

B=nd.*(1+nd)+(sigma.*c.^2.*(-1-2.*m.^2+2.*nd.*(1+nd)))./(-3+4.*nd.*(1+nd))-(ev+sigma.*c.^2);

% C=(c.^2.*(m-nd).*(1+m-nd))./(3+4.*(-2+nd).*nd);
C=sigma.*c.^2.*sqrt(((m-nd).*(1+m-nd).*(-1+m+nd).*(m+nd))./((1-2.*nd).^2.*(-3+2.*nd).*(1+2.*nd)));
% C=sqrt((c.^4.*(m-nd).^3.*(1+m-nd).^3.*(1+2.*nd))./((1-2.*nd).^2.*(-1+m+nd).*(m+nd).*(-3+2.*nd).^3));

% nd=nd+2;
% A=(c.^2.*(1+m+nd).*(2+m+nd))./((3+2.*nd).*(5+2.*nd))
A=sigma.*c.^2.*sqrt(((-2+m-nd).*(-1+m-nd).*(1+m+nd).*(2+m+nd))./((1+2.*nd).*(3+2.*nd).^2.*(5+2.*nd)));
% A=sqrt((c.^4.*(1+m+nd).^3.*(2+m+nd).^3.*(1+2.*nd))./((-2+m-nd).*(-1+m-nd).*(3+2.*nd).^2.*(5+2.*nd).^3));

% C1=sqrt(C.*A);
u1u0=C./(-B-A.*u2u1);
end

function u1u2=up_recursion_ratio_norm(sigma,nd,m,c,ev,u0u1);
% note that this is u_nd./u_nd+2. input is u_nd-2./u_nd
%note: no handling of the zero.

B=nd.*(1+nd)+(sigma.*c.^2.*(-1-2.*m.^2+2.*nd.*(1+nd)))./(-3+4.*nd.*(1+nd))-(ev+sigma.*c.^2);

C=(sigma.*c.^2.*(m-nd).*(1+m-nd))./(3+4.*(-2+nd).*nd)./sqrt(((m-nd).*(1+m-nd).*(1+2.*nd))./((-1+m+nd).*(m+nd).*(-3+2.*nd)));
% C=(c.^2).*sqrt(((m-nd).*(1+m-nd).*(-1+m+nd).*(m+nd))./((-3+2.*nd).*(1+2.*nd)))./(-1+2.*nd);

A=(sigma.*c.^2.*(1+m+nd).*(2+m+nd))./((3+2.*nd).*(5+2.*nd)).*sqrt(((-2+m-nd).*(-1+m-nd).*(5+2.*nd))./((1+m+nd).*(2+m+nd).*(1+2.*nd)));
% A=(c.^2).*sqrt(((-2+m-nd).*(-1+m-nd).*(1+m+nd).*(2+m+nd))./((1+2.*nd).*(5+2.*nd)))./(3+2.*nd);

u1u2=-A./(B+C.*u0u1);

end
