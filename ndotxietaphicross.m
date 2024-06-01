function [J_element]=ndotxietaphicross(n,A,B,ds);
% ndotxietaphicross - computes the J_element. A helper function to make
% 						sure that no thought be needed.
%
% usage:
%
% [J_element]=ndotxietaphicross(n,A,B,ds);
%
% J_element -- the result of sum(ds*n.(A x B)). Hopefully this is an integral.
% n         -- surface normal for each point.
% A         -- numpoints x 3-vector.
% B         -- numpoints x 3-vector.
% dS        -- surface element, or other scalar factor needed for each point.
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

J_element=sum(ds.*(n(:,1).*(A2B3-A3B2)+n(:,2).*(A3B1-A1B3)+n(:,3).*(A1B2-A2B1)));
