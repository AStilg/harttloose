function [yn,dyn] = sbessely(n,kr)
% sbesselh1 - spherical hankel function hn(kr) of the first kind
%
% Usage:
% hn = sbesselh1(n,kr)
%
% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
%
% See besselh for more details
%
% PACKAGE INFO

kr=kr(:);
n=n(:);

if nargout==2
    n=[n;n+1];
end

[n,kr]=meshgrid(n,kr);
[yn] = bessely(n+1/2,kr);

small_args = find( abs(kr) < 1e-15 );
not_small_args = find( ~(abs(kr) < 1e-15) );

if length(kr) == 1 && abs(kr) < ((1e-15))

    yn = -kr.^(-n-1) .*repmat(prod(1:2:(2*n+1)),[1,length(n)]);
    
elseif length(kr) == 1 && ~((1e-15).^(1./(-n-1)))

    yn = sqrt(pi./(2*(kr))) .* yn;
elseif length(n) == 1

    yn(not_small_args) = ...
        sqrt(pi./(2*(kr(not_small_args)))) .* yn(not_small_args);
    yn(small_args) = -kr(small_args).^(-n-1) .* prod(1:2:(2*n+1));
else % both n and kr are vectors

    yn(not_small_args) = ...
        sqrt(pi./(2*(kr(not_small_args)))) .* yn(not_small_args);
    yn(small_args) = -kr(small_args).^(-n(small_args)-1) .* ...
        prod(1:2:(2*n(small_args)+1));
end

yn=(yn);

if nargout==2
    dyn=((n(1:end,1:end/2)+1) .* yn(1:end,1:end/2)./kr(1:end,1:end/2)-yn(1:end,end/2+1:end));
    yn=yn(1:end,1:end/2);
end



return
