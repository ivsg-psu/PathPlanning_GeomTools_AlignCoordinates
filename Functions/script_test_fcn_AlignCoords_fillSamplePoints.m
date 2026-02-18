% script_test_fcn_AlignCoords_fillSamplePoints
% Tests function: fcn_AlignCoords_fillSamplePoints

% 
% REVISION HISTORY:
% 
% 2023_03_23 by Sean Brennan
% - first write of function


close all;



%% Try basic call 
% An empty point type gives default
figNum = 1;
fcn_AlignCoords_fillSamplePoints([], figNum);
assert(true);

%% Show all the options
% Fills in some sample points, in homogenous form
figNum = 2;
figure(figNum);
clf;

for type_of_points = 1:4
    subplot(2,2,type_of_points);
    start_points = fcn_AlignCoords_fillSamplePoints(type_of_points, figNum);
end
sgtitle('Demonstration of fcn_AlignCoords_generate2DTransformMatrix', 'Interpreter', 'none');

%% Fail conditions
if 1==0
    %% Bad input - wrong number of arguments
    fcn_AlignCoords_fillSamplePoints(24,25);
    
end
