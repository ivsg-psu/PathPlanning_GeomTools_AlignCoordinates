% script_test_fcn_AlignCoords_regressionFitScaleFactor
% Tests function: fcn_AlignCoords_regressionFitScaleFactor

% 
% REVISION HISTORY:
% 
% 2023_03_23 by Sean Brennan
% -- first write of function


close all;

% Fill in a sample transform matrix parameters
S = 2;
tx = 2;
ty = 7;
theta = 20*pi/180;


% Fill in transformation matrix
order_string = 'rts';
T = fcn_AlignCoords_generate2DTransformMatrix( ...
    S, theta, tx, ty, order_string);

% Fills in some sample points, in homogenous form
start_points = fcn_AlignCoords_fillSamplePoints;

% apply the transform
moved_points = (T*start_points')';

coord_base_points = start_points;
coord_xform_points = moved_points;


%% Try basic call
fig_num = 1;
S_calculated = fcn_AlignCoords_regressionFitScaleFactor(coord_base_points(:,1:2), coord_xform_points(:,1:2), fig_num);
assert(abs(S-S_calculated)<1E-10);


%% Fail conditions
if 1==0
    %% Bad input - wrong number of arguments
    fcn_AlignCoords_generate2DTransformMatrix(24);
    
    URHERE
    %% Bad input - wrong number of arguments
    fcn_CodeX_generateNumbersLike_KEEP(24,3,4);
    %% Bad input - not integer
    fcn_CodeX_generateNumbersLike_KEEP(2.34,5);
    %% Bad input - not positive
    fcn_CodeX_generateNumbersLike_KEEP(-2,5);
    %% Bad input - not integer
    fcn_CodeX_generateNumbersLike_KEEP(2,5.4);
    %% Bad input - not positive
    fcn_CodeX_generateNumbersLike_KEEP(2,-5.4)
end
