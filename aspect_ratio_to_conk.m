function [isProlate,conk]=aspect_ratio_to_conk(ac);
% aspect_ratio_to_conk - takes the z-axis and x-axis maximum object size
% 						and estimates a "good" interfocal distance. This
%                       is a helper function.
% 
% PACKAGE INFO

if numel(ac)==2
    ac=ac(:);
end

major_axis=ac(2,:);
minor_axis=ac(1,:);
aspect_ratio=major_axis./minor_axis;

isProlate=zeros(1,size(ac,2));
isProlate(:)=true;
isProlate(ac(2,:)<ac(1,:))=false;

sigma=2*isProlate-1;

conk=(sqrt(sigma.*(aspect_ratio.^2-1)).*minor_axis); % conk_m