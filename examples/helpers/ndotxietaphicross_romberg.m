function [J_element,to_rearrange]=ndotxietaphicross_romberg(n,A,B,all_coordinate,surfs);
% ndotxietaphicross_romberg - computes the J_element. A helper function to make
% 						sure that no thought be needed. Note: this was developed 
%                       for a romberg integration.
%
% usage:
%
% [J_element,to_rearrange]=ndotxietaphicross_romberg(n,A,B,all_coordinate,surfs);
%
% J_element -- the result of romberg integration. Hopefully this is an integral.
% n         -- surface normal for each point.
% A         -- numpoints x 3-vector.
% B         -- numpoints x 3-vector.
% all_coordinate -- For every connected surface concat the coordinates. Each surface
%                   must have EQUAL power of 2 points.
% surfs     -- numer of surfaces in the integration.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

A2B1=A(:,2).*B(:,1);
A3B1=A(:,3).*B(:,1);
A1B2=A(:,1).*B(:,2);
A3B2=A(:,3).*B(:,2);
A1B3=A(:,1).*B(:,3);
A2B3=A(:,2).*B(:,3);

to_rearrange=(n(:,1).*(A2B3-A3B2)+n(:,2).*(A3B1-A1B3)+n(:,3).*(A1B2-A2B1));
to_romberg=reshape(to_rearrange,[],surfs);
to_coordinates=reshape(all_coordinate,[],surfs);
J_element=sum(romberg_int_grid(to_romberg,to_coordinates));
