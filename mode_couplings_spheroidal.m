function modeMatrix=mode_couplings_spheroidal(Nmax,symmetries);
%mode_couplings.m : Construct a matrix which predicts the non-zero sT-matrix
%                   elements for particular symmetries.
%
% USAGE:
%
% smodeMatrix=mode_couplings_spheroidal(Nmax,symmetries);
%
% symmetries = [m_coupling,xymirror]
% m_couling = minimum difference between m and m'
% xymirror  = does the particle have mirror symmetry in xy-plane
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

% %%%%%%%%%%TEST CASES%%%%%%%%%%
% Nmax=5;
% symmetries=[5,1];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isinf(symmetries(1))
    symmetries(1)=2*Nmax+1;
end

rotSym=symmetries(1);
mirrorSym=symmetries(2);

[n,m]=combined_index([0:Nmax*(Nmax+2)].');

modeMatrix=spalloc(Nmax*(Nmax+2)*2+2,Nmax*(Nmax+2)*2+2,((Nmax*(Nmax+2)+1)*2)^2/2);

if mirrorSym
    Z=(rem(n+m,2)==0);
else
    Z=ones(size(n));
end
Z2=Z;

%rotational couplings
for ii=1:length(n)
    %indexes this column

    if mirrorSym;
        if rem(n(ii)+ii-1,2)
            Z2=~Z;
        else
            Z2=Z;
        end
    end
    
    X=all([Z2,~rem(rotSym-m+m(ii),rotSym)],2);
    
    if mirrorSym
        Y=all([~Z2,~rem(rotSym-m+m(ii),rotSym)],2);
    else
        Y=X;
    end

    if m(ii)==0
        Y(m==0)=0;
    end
    
    cix=n(X).*n(X)+n(X)+m(X)+1;
    ciy=n(Y).*n(Y)+n(Y)+m(Y)+1;
    
    modeMatrix(cix,ii)=1;
    modeMatrix(Nmax*(Nmax+2)+1+cix,Nmax*(Nmax+2)+1+ii)=1;
    modeMatrix(ciy,Nmax*(Nmax+2)+1+ii)=1;
    modeMatrix(Nmax*(Nmax+2)+1+ciy,ii)=1;
    
%     if mirrorSym%
%         Z=~Z;
%     end
    
end
