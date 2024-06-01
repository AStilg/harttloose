function [hn,dhn] = sbesselh1(n,kr)
% sbesselh1 - spherical hankel function hn(kr) of the first kind
%
% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
%
% Usage:
%
% hn = sbesselh1(n,z);
% OR
% [hn,dzhn] = sbesselh1(n,z);
%
% where dzhn is the derivative of the appropriate Ricatti-Bessel function
% divided by z.
%
% See besselj and bessely for more details
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

% if round(kr)>=n
%     if nargout==2
%         n=[n;n+1];
%     end
%     
%     [n,kr]=meshgrid(n,kr);
%     [hn] = besselh(n+1/2,1,kr);
%     
%     hn = sqrt(pi./(2*kr)) .* (hn);
%     
%     if nargout==2
%         dhn=((n(1:end,1:end/2)+1) .* hn(1:end,1:end/2)./kr(1:end,1:end/2)-hn(1:end,end/2+1:end));
%         hn=hn(1:end,1:end/2);
%     end
% else
    [jn,djn]=sbesselj(n,kr);
    [yn,dyn]=sbessely(n,kr);
    hn=jn+1i*yn;
    if nargout==2
        dhn=djn+1i*dyn;
    end
% end

return
