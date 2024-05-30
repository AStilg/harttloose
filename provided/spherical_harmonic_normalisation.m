function [N,logflag]=spherical_harmonic_normalisation(n,m);
% spherical_harmonic_normalisation - computes the normalisation factor, N,
% 									such that \bar{P}_{nm}=N_{nm}*P_{nm}. 
%									Y_{nm}=\bar{P}_{nm} e^{im\phi}.
%
% USAGE:
%
% [N,logflag]=spherical_harmonic_normalisation(n,m)
%
% n -- radial mode
% m -- azimuthal mode
%
% PACKAGE INFO

logflag=0;

if length(n)<length(m)
    n=ones(size(m))*n;
end
if length(m)<length(n)
    m=ones(size(n))*m;
end

N=0*n;
for ii=1:length(n)
    
    nmm=[1:n(ii)-m(ii)];
    npm=[1:n(ii)+m(ii)];
    
    mindex=min(length(nmm),length(npm));
    
    common_product=prod(sqrt(nmm(1:mindex)./npm(1:mindex)));
    if m(ii)+n(ii)>n(ii)-m(ii)
        common_product=common_product.*prod(1./sqrt(npm(mindex+1:end)));
    end
    if n(ii)-m(ii)>n(ii)+m(ii)
        common_product=common_product.*prod(sqrt(nmm(mindex+1:end)));
    end
    
    N(ii)=(-1).^m(ii).*sqrt((2*n(ii)+1)/(4*pi)).*common_product;
    if m(ii)==0
        N(ii)=(-1).^m(ii).*sqrt((2*n(ii)+1)/(4*pi));
    end
    
    if or(any(N(ii)==0),any(isnan(N(ii))))
        logflag=1;
        warning('At least one number too small/big... going to logs')
        N(ii)=log10(sqrt((2*n+1)/(4*pi)))+.5*(logfactorial(n(ii)-m(ii))-logfactorial(n(ii)+m(ii)));
    end
end
