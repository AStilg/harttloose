function [jn,djn,d2jn] = sbesselj(n,kr)
% sbesselj - spherical bessel function jn(kr)
%
% jn(kr) = sqrt(pi/2kr) Jn+0.5(kr)
%
% Usage:
%
% jn = sbessel(n,z);
% OR
% [jn,dzjn] = sbessel(n,z);
%
% where dzjn is the derivative of the appropriate Ricatti-Bessel function
% divided by z.
%
% See besselj for more details
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

kr=kr(:);
n=n(:);

if nargout>1
    n=[n;n+1];
end

[n,kr]=meshgrid(n,kr);
% [jn] = (sign(kr)).^(n).*besselj(n+1/2,abs(kr));
[jn] = besselj(n+1/2,(kr));

small_args = find( abs(kr) < 1e-15 );
not_small_args = find( ~(abs(kr) < 1e-15) );

if length(kr) == 1 && abs(kr) < 1e-15

    jn = kr.^n ./repmat(prod(1:2:(2*n+1)),[1,length(n)]);
    
elseif length(kr) == 1 && ~(abs(kr) < 1e-15)

    jn = sqrt(pi./(2*(kr))) .* jn;
elseif length(n) == 1

    jn(not_small_args) = ...
        sqrt(pi./(2*(kr(not_small_args)))) .* jn(not_small_args);
    jn(small_args) = kr(small_args).^n ./ prod(1:2:(2*n+1));
else % both n and kr are vectors

    jn(not_small_args) = ...
        sqrt(pi./(2*(kr(not_small_args)))) .* jn(not_small_args);
    jn(small_args) = kr(small_args).^n(small_args) ./ ...
        prod(1:2:(2*n(small_args)+1));
end
jn=(jn);

if nargout>1
    small_args = find( abs(kr(1:end,1:end/2)) < 1e-15 );
    djn=((n(1:end,1:end/2)+1) .* jn(1:end,1:end/2)./kr(1:end,1:end/2)-jn(1:end,end/2+1:end));
    d2jn=(n(1:end,1:end/2).*(n(1:end,1:end/2)-1)./kr(1:end,1:end/2).^2-1).*jn(1:end,1:end/2)+2.* jn(1:end,end/2+1:end)./kr(1:end,1:end/2);
    if length(small_args)>0
        djn(small_args)=n(small_args).*kr(small_args).^(n(small_args)-1) ./ prod(1:2:(2*n(small_args)+1));
        nt=n(small_args);
        nt(nt<2)=2;
        d2jn(small_args)=(n(small_args).*(n(small_args)-1).*kr(small_args).^(nt-2)+(n(small_args)+2).*(n(small_args)-1).*kr(small_args).^(n(small_args))./(2*n(small_args)+3)) ./ prod(1:2:(2*n(small_args)+1));
    end
    djn(isnan(djn))=0;
    d2jn(isnan(d2jn))=0;
    jn=jn(1:end,1:end/2);
end

return
