function [I,R]=romberg_int_grid(f,x); 
% romber_int_grid - calculates romberg integration to however many
%                   nested Richardson extrapolations it can.
%
% USAGE:
% [I,R]=romberg_int_grid(f,x)
% 
% x  -- coordinate points, must be length=2^n+1;
% f  -- function points
%
% f, and, x may be equal size matrix of [n-points x integrals]
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if nargin==1
    x=([1:size(f,1)]').*ones(1,size(f,2));
end

N=(log(size(x,1)-1)/log(2));

if length(x)~=2^floor(N)+1
    error('length of data must be 2^N+1 only');
    % return;
end

N=round(N);

R=zeros(N+1,size(f,2),N+1);

R(1,:,1)=(f(1,:)+f(end,:)).*(x(end,:)-x(1,:))/2;

%% trapezoidal step
for ii=1:N
    g=f(1:2^N/2^ii:end,:);
    R(ii+1,:,1)=sum(diff(x(1:2^N/2^ii:end,:),1,1) .* (g(1:end-1,:) + g(2:end,:))/2,1);
end

%%
%romberg richardson extrapolation
for jj=2:N+1
    R(jj:end,:,jj)=((4^(jj-1))*R(jj:end,:,jj-1)-R(jj-1:end-1,:,jj-1))/(4^(jj-1)-1);
end
I=R(end,:,end);
