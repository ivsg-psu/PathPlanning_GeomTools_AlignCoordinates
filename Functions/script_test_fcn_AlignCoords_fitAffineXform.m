% script_test_fcn_AlignCoords_fitAffineXform
% Tests function: fcn_AlignCoords_fitAffineXform

% 
% REVISION HISTORY:
% 
% 2023_03_23 by Sean Brennan
% -- first write of function


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

fig_num = 1;
[T_calculated,err] = fcn_AlignCoords_fitAffineXform(coord_base_points(:,1:2), coord_xform_points(:,1:2), fig_num); % Find optimal transform

% Is the error small?
assert(max(err,[],'all')<1E-10);

% Show that T_calculated will not be the same as T
fprintf(1,'The transformation calculated via affine will not be the same as T. \n');
fprintf(1,'To illustrate, here is the result with no noise: \n');
fprintf(1,'This is T used to perform coordinate rotations:\n');
disp(T)
fprintf(1,'This is T_calculated resulting from the affine fit:\n');
disp(T_calculated)
fprintf(1,'Note: the fit is very good as the errors are very small, with maximum error: %s\n',num2str(max(err,[],'all'),6));




%% Test with noisy data
fig_num = 2;
figure(fig_num);
clf;
hold on;

% add noise now to the moved points
coord_base_points = start_points(:,1:2);
offsets = 0.5*randn(length(start_points(:,1)),2);
coord_xform_points = moved_points(:,1:2) + offsets;
offset_distances = sum(offsets.^2,2).^0.5;

% Calculate the result
[T_calculated,err] = fcn_AlignCoords_fitAffineXform(coord_base_points(:,1:2), coord_xform_points(:,1:2), fig_num); % Find optimal transform

fprintf(1,'Now, here is the affine result with some noise: \n');
fprintf(1,'This is T used to perform coordinate rotations:\n');
disp(T)
fprintf(1,'This is T_calculated resulting from the affine fit:\n');
disp(T_calculated)
fprintf(1,'Note: the fit is poor even though the errors are small, with maximum error: %s\n',num2str(max(err,[],'all'),6));


[T_calculated,~,~,~,err] = fcn_AlignCoords_fit2DCoordinates(coord_base_points(:,1:2), coord_xform_points(:,1:2)); 
sgtitle('Demonstration of fcn_AlignCoords_fit2DCoordinates', 'Interpreter', 'none','FontSize',12);

fprintf(1,'Finally, here is the 2D coordinate fit result with the same inputs as affine: \n');
fprintf(1,'This is T used to perform coordinate rotations:\n');
disp(T)
fprintf(1,'This is T_calculated resulting from the 2D rotation fit:\n');
disp(T_calculated)
fprintf(1,'Note: the fit is still poor, the rotation is preserved, however. It has maximum error: %s\n',num2str(max(err,[],'all'),6));



fprintf(1,'Errors due to scaling:\n');
disp(offset_distances - err);

   
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
