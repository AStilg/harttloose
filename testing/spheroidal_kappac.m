function [knm1]=spheroidal_kappac(isProlate,n,m,c,allvalues);
% spheroidal_kappac - calculates the ratio between radial and angular spheroidal
%                     wavefunctions for specified c. Only for S1 and R1.
%
% usage :
%
% knmc=spheroidal_kappac(isProlate,n,m,c,allvalues);
%
% isProlate - true/false
% n         - pseudo radial coordinate
% m         - pseudo angular coordinate
% c         - interfocal distance
% allvalues - true/false (all modes computed of same parity)
% knmc      - S1nmc/R1nmc
% 
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if nargin==4
    allvalues=false;
end

if isProlate
    knm1=(spheroidalS1(isProlate,n,abs(m),c,sqrt(2)*1i,allvalues)./spheroidalR1(isProlate,n,abs(m),c,sqrt(2)*1i,allvalues));
    if mod(m,2)
    	knm1=-1i*knm1;
    end
else
    knm1=(spheroidalS1(isProlate,n,abs(m),c,[1.01i,-1e-8i],allvalues)./spheroidalR1(isProlate,n,abs(m),c,[-1.01i,1e-8]',allvalues));
    if mod(m,2)
    	knm1=-1i*knm1;
    end
    nd=n;
    if size(knm1,1)~=1
        nd=[0:size(knm1,1)-1]'*2+abs(m)+mod(n+m,2);
    end
    [~,indx]=max(abs(knm1),[],2); %maximum value is going to be correct.
    switchcol=(indx-1)*size(knm1,1)+[1:size(knm1,1)]';
    knm1=knm1(switchcol);
end
knm1=1./knm1;
