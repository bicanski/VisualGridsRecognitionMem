

% Code developed by Andrej Bicanski 
% andrej.bicanski@gmail.com
% www.andrejbicanski.com
%
% model published in Current Biology
%
% Bicanski A, Burgess N. - A computational model of recognition 
% memory via grid cells. Current Biology, 2019, 29, 1–12. 
% DOI: 10.1016/j.cub.2019.01.077
%
% This script calculates the x and y distance given start and goal grid
% cell population vectors, please see Bush et al. 2015 for details on the
% distance cell model. Note, any equivalent vector computation could be
% be substituted for the distance model without changes to the rest of the model.
%
% Bush D, Barry C, Manson D, Burgess N. Using grid cells for navigation. 
% Neuron. 2015 Aug 5;87(3):507-20.



function [x_diff,y_diff] = DOE_subr_GCs2Dist(cGCs,tGCs,AllPosGC2Dy_wts,AllPosGC2Dx_wts, ... 
                                    yDgoal2yUPrdout_wts,yDgoal2yDOrdout_wts,xDgoal2yUPrdout_wts,xDgoal2yDOrdout_wts, ...
                                           yDcurr2yUPrdout_wts,yDcurr2yDOrdout_wts,xDcurr2yUPrdout_wts,xDcurr2yDOrdout_wts);

Dy_curr_rates = AllPosGC2Dy_wts * cGCs;                % look up distance rates at target and current location
Dy_goal_rates = AllPosGC2Dy_wts * tGCs;
Dx_curr_rates = AllPosGC2Dx_wts * cGCs;
Dx_goal_rates = AllPosGC2Dx_wts * tGCs;

Dy_curr_rates(Dy_curr_rates<max(Dy_curr_rates)) = 0;   % simple winner take all
Dy_goal_rates(Dy_goal_rates<max(Dy_goal_rates)) = 0;
Dx_curr_rates(Dx_curr_rates<max(Dx_curr_rates)) = 0;
Dx_goal_rates(Dx_goal_rates<max(Dx_goal_rates)) = 0;

Dy_curr_rates = min(1,Dy_curr_rates);%/max(Dy_curr_rates);      % max 1
Dy_goal_rates = min(1,Dy_goal_rates);%/max(Dy_goal_rates);
Dx_curr_rates = min(1,Dx_curr_rates);%/max(Dx_curr_rates);
Dx_goal_rates = min(1,Dx_goal_rates);%/max(Dx_goal_rates);

yUPrdout_rate = yDcurr2yUPrdout_wts * Dy_curr_rates + yDgoal2yUPrdout_wts * Dy_goal_rates;    % readout
yDOrdout_rate = yDcurr2yDOrdout_wts * Dy_curr_rates + yDgoal2yDOrdout_wts * Dy_goal_rates;
xUPrdout_rate = xDcurr2yUPrdout_wts * Dx_curr_rates + xDgoal2yUPrdout_wts * Dx_goal_rates;
xDOrdout_rate = xDcurr2yDOrdout_wts * Dx_curr_rates + xDgoal2yDOrdout_wts * Dx_goal_rates;

x_diff = (xUPrdout_rate-xDOrdout_rate)/2;   % map readout rate to position on plane
y_diff = (yUPrdout_rate-yDOrdout_rate)/2;

if x_diff==0 && y_diff==0
    test=1;
end

