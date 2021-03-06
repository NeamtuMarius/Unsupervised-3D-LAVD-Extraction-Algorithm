%% References:
% [1]: G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn,
%      Defining coherent vortces objectively from the vorticity. submitted (2015)
% [2]: G. Haller, Dynamically consistent rotation and stretch tensors for
%      finite continuum deformation. submitted (2015).
%%
% function flag = IsContourConvex(xt,yt,DeficiencyThresh)

% Input arguments:
%   xt: x-component of a contour
%   yt: y-component of a contour
%   DeficiencyThresh: convexity deficiency threshold - as a percentage of
%   convex hull

% Output arguments:
%   flag: determines if a contour is convex - true
%--------------------------------------------------------------------------
% Author: Alireza Hadjighasem  alirezah@ethz.ch
% http://www.zfm.ethz.ch/~hadjighasem/index.html
%--------------------------------------------------------------------------
function flag = IsContourConvex(xt,yt,DeficiencyThresh)
pointsAreCollinear = @(xy) rank(xy(2:end,:) - xy(1,:)) == 1; % function to check collinearity in 2d
if length(xt) < 3 | pointsAreCollinear([xt yt]) == 1 % if points are collienar or only only 2
    flag = false;
    return
end

[~,Hull] = convhull(xt,yt);
Area = polyarea(xt,yt);
AreaDeficiency = abs(Hull-Area)./Area*100;

if AreaDeficiency <= DeficiencyThresh
    flag = true;
else
    flag = false;
end


end



