figure; % Set up a new figure
hold on; % Hold the plots

% Define a transformation
S = 2;
theta = 45*pi/180;
tx = 3;
ty = 7;
order_string = 'tsr';
T = fcn_AlignCoords_generate2DTransformMatrix(S, theta, tx, ty, order_string);

Npoints = 10; % How many points should we demo with?
coord_base_points = rand(Npoints,2); % Set up some test points
plot(coord_base_points(:,1),coord_base_points(:,2),'b.','Markersize',10); % Plot the result

% Convert to homogenous form in prep for transform
colOfOne = ones(Npoints,1); % Define column of 1's
normalized_coord_base_points = [coord_base_points(:,1:2) colOfOne];  % Make homogenous
normalized_coord_xform_points = (T*normalized_coord_base_points')';

% Plot the results
plot(normalized_coord_xform_points(:,1),normalized_coord_xform_points(:,2),'r.','Markersize',10); % Plot the result

% Invert them back
moved_back_points = (T\normalized_coord_xform_points')';

% Plot the results, and show that the moved back points, given by green
% circles, land back at the original points, shown in blue
plot(moved_back_points(:,1),moved_back_points(:,2),'go','Markersize',10); % Plot the result