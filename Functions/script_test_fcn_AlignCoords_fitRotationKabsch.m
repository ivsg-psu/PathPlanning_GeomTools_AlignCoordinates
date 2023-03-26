% script_test_fcn_AlignCoords_fitRotationKabsch
% Tests function: fcn_AlignCoords_fitRotationKabsch

% 
% REVISION HISTORY:
% 
% 2023_03_23 by Sean Brennan
% -- first write of function


close all;

% Fill in a sample transform matrix parameters
S = 1;
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
[R_calculated,t_rotated,err] = fcn_AlignCoords_fitRotationKabsch(coord_base_points(:,1:2), coord_xform_points(:,1:2), fig_num); % Find optimal transform


% Calculate unrotated version of t, and in column form. Note: it's in the
% negative direction because the fitting method is from moved to base.
t_calculated = -(t_rotated/R_calculated)';

% Is the error small?
assert(max(err,[],'all')<1E-10,'Error did not pass assertion');

% Does the rotation match?
assert(max(R_calculated-T(1:2,1:2),[],'all')<1E-10,'Rotation did not pass assertion');

% Does the translation match?
assert(max(t_calculated-T(1:2,3),[],'all')<1E-10,'Translation did not pass assertion');


%% Test across dimensions
fig_num = 2;
figure(fig_num);
clf;
hold on;

% Number of points
Numpoints = 10; %randi([1 1e2]); % Number of points (integer > 0)

% m is the dimension
for m = 1:4
    subplot(2,2,m);
    
    R_orig = orth(randn(m));                % Random rotation matrix
    t_orig = randn(1,m);                    % Random translation vector
    
    coord_base_points = rand(Numpoints,m);                 % Random point cloud
    coord_xform_points = bsxfun(@plus,coord_base_points*R_orig,t_orig);   % Transform points

    [R_calculated,t_rotated,err] = fcn_AlignCoords_fitRotationKabsch(coord_base_points, coord_xform_points); % Find optimal transform
    
    P1 = coord_xform_points;
    P2 = coord_base_points;
    P3 = bsxfun(@plus,P1*R_calculated,t_rotated);              % Apply optimal transform to points
    maxerr = max(err);                      % Maximum RMSE
    
    % Plot points
    m = length(P1(1,:)); % What is the dimension of the points?
    n = length(P1(:,1)); % How many points?
    
    
    
    % Change the plotting depending on dimension
    if m == 1
        z = zeros(n,1);
        plot(1+z,P1(:,1),'ro',...
            2+z,P2(:,1),'bo',...
            2+z,P3(:,1),'k.');
    elseif m == 2
        plot(P1(:,1),P1(:,2),'ro',...
            P2(:,1),P2(:,2),'bo',...
            P3(:,1),P3(:,2),'k.');
    else
        plot3(P1(:,1),P1(:,2),P1(:,3),'ro',...
            P2(:,1),P2(:,2),P2(:,3),'bo',...
            P3(:,1),P3(:,2),P3(:,3),'k.');
        view(3);
    end
    
    
    axis equal;
    grid on;
    legend('Moved points','Target points','Moved points Transformed');
    title([int2str(m) ' dimensions, ' int2str(n) ' points, Maximum RMSE: ' ...
        num2str(maxerr,6)]);
    
    
    
    t_calculated = -t_rotated/R_calculated;
    
    % Is the error small?    
    assert(max(err,[],'all')<1E-10);

    % Does the rotation match? (Note: rotations in 3D are non-unique, so we
    % can check the eigenvalues)
    diff = sort(eig(R_orig),'ComparisonMethod','real')-sort(eig(R_calculated),'ComparisonMethod','real');
    assert(max(abs(diff),[],'all')<1E-10);
    
    % Does the translation match?
    assert(max(t_calculated-t_orig,[],'all')<1E-10);

end

%% Test with noisy data
fig_num = 3;
figure(fig_num);
clf;

hold on;

% Fills in some sample points, in homogenous form
start_points = fcn_AlignCoords_fillSamplePoints;

% Fill in a sample transform matrix parameters
S = 1;
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
coord_xform_points = moved_points(:,1:2) + 0.1*randn(length(start_points(:,1)),2);

[R_calculated,t_rotated,err] = fcn_AlignCoords_fitRotationKabsch(coord_base_points, coord_xform_points); % Find optimal transform

% Calculate unrotated version of t, and in column form. Note: it's in the
% negative direction because the fitting method is from moved to base.
t_calculated = -t_rotated/R_calculated;

P3 = bsxfun(@plus,coord_xform_points*R_calculated,t_rotated);
plot(coord_xform_points(:,1),coord_xform_points(:,2),'r.-','LineWidth',3,'MarkerSize',20);
plot(coord_base_points(:,1),coord_base_points(:,2),'b.-','LineWidth',3,'MarkerSize',20);
plot(P3(:,1),P3(:,2),'ko-','LineWidth',3,'MarkerSize',10);

axis equal;
grid on;
legend('Moved points','Target points','Moved points Transformed');
m = 2;
Npoints = length(coord_base_points(:,1));
title([int2str(m) ' dimensions, ' int2str(Npoints) ' points, Maximum RMSE: ' ...
    num2str(max(err),6)]);



   
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
