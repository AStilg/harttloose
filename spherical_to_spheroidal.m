function A_mod = spherical_to_spheroidal(isProlate,c,A,nmaxout,force_vector);
% spherical_to_spheroidal - transforms TE or TM coefficients to spheroidal
%							equivalents. This is helper function.
%
% usage:
%
% A_mod = spherical_to_spheroidal(isProlate,c,A);
% A_mod = spherical_to_spheroidal(isProlate,c,A,nmaxout);
%
% where:
%
% isProlate -- true/false
% c         -- interfocal number. Note: two mediums of operation can be set
%              in the order of [c_out,c_in].
% A         -- vector or matrix of either TE or TM coefficients.
% nmaxout   -- makes the output nmax different to input.
%
% Note 1: when the A is a matrix it is a block of either TETE, TETM, TMTE, 
% or TMTM  couplings. If a whole T-matrix is entered the result will be invalid.
%
% Note 2: this is a FUNDAMENTALLY LOSSY process. You will need at least an
% additional nmax radial index to preserve information for multiple
% transformations.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

nmax=floor(sqrt(size(A,1)));
if nargin<4
	nmaxout=nmax;
end
if nargin<5
	force_vector=false;
end

total_index_in=[1:nmax*(nmax+2)]+1;
total_index_out=[1:(nmaxout)*((nmaxout)+2)+1];

[U,Nnm0,~,im_v]=spheroidal_expansion(isProlate,c(end),(max([nmax,nmaxout])));

Umod=U(total_index_out,total_index_in).*(im_v(total_index_out).*im_v(total_index_in)');
N=Nnm0(total_index_in);

if numel(c)>1
    [U,~,~,im_v]=spheroidal_expansion(isProlate,c(1),(max([nmax,nmaxout])));
    Umod2=U(total_index_out,total_index_in).*(im_v(total_index_out).*im_v(total_index_in)');
    
    if and(size(A,2)>1,~force_vector)
        A_mod=Umod2*((A.*N.'./N)*Umod');
    else
        error('A vector cannot have two mediums of operation!');
    end    
else
    
    if and(size(A,2)>1,~force_vector)
        A_mod=Umod*((A.*N.'./N)*Umod');
    else
        A_mod=(Umod)*(A./N);
    end
    
end
