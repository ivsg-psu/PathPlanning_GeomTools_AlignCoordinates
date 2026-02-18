% script_test_fcn_AlignCoords_fit2DCoordinates
% Tests function: fcn_AlignCoords_fit2DCoordinates

% 
% REVISION HISTORY:
% 
% 2023_03_23 by Sean Brennan
% - first write of function


close all;

% Fill in a sample transform matrix parameters
S = 2;
tx = 10;
ty = 15;
theta = 50*pi/180;


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

figNum = 1;
[T_calculated,R_calculated,S_calculated,t_calculated,err] = fcn_AlignCoords_fit2DCoordinates(coord_base_points(:,1:2), coord_xform_points(:,1:2), figNum); % Find optimal transform

% Is the error small?
assert(max(err,[],'all')<1E-10);

% Does the T matrix match?
assert(max(T_calculated-T,[],'all')<1E-10, 'Transformation matrix, T, did not match');

% Does the scaling match?
assert(max(S_calculated-S,[],'all')<1E-10, 'Scaling, S, did not match');

% Does the rotation match?
assert(max(R_calculated-T(1:2,1:2)/S,[],'all')<1E-10, 'Rotation matrix, R, did not match');

% Does the translation match?
assert(max(t_calculated-T(1:2,3)/S,[],'all')<1E-10, 'translation did not match');




%% Test with noisy data
figNum = 2;
figure(figNum);
clf;

hold on;

% Fills in some sample points, in homogenous form
start_points = fcn_AlignCoords_fillSamplePoints;

% Fill in a sample transform matrix parameters
S = 3.2;
tx = 2;
ty = 7;
theta = 50*pi/180;


% Fill in transformation matrix
order_string = 'rts';
T = fcn_AlignCoords_generate2DTransformMatrix( ...
    S, theta, tx, ty, order_string);

% apply the transform, BUT WITH NOISE
moved_points = (T*start_points')';

coord_base_points = start_points(:,1:2);
offsets = 0.5*randn(length(start_points(:,1)),2);
coord_xform_points = moved_points(:,1:2) + offsets;
offset_distances = sum(offsets.^2,2).^0.5;

% Calculate the result
[T_calculated,R_calculated,S_calculated,t_calculated,err] = fcn_AlignCoords_fit2DCoordinates(coord_base_points(:,1:2), coord_xform_points(:,1:2), figNum); % Find optimal transform

fprintf(1,'T true: \n');
disp(T);

fprintf(1,'T calculated: \n');
disp(T_calculated);

fprintf(1,'Errors due to scaling:\n');
disp(offset_distances - err);


% % Is the error small?
% assert(max(err,[],'all')<1E-10);
% 
% % Does the T matrix match?
% assert(max(T_calculated-T,[],'all')<1E-10, 'Transformation matrix, T, did not match');
% 
% % Does the scaling match?
% assert(max(S_calculated-S,[],'all')<1E-10, 'Scaling, S, did not match');
% 
% % Does the rotation match?
% assert(max(R_calculated-T(1:2,1:2)/S,[],'all')<1E-10, 'Rotation matrix, R, did not match');
% 
% % Does the translation match?
% assert(max(t_calculated-T(1:2,3)/S,[],'all')<1E-10, 'translation did not match');
% 



   
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
