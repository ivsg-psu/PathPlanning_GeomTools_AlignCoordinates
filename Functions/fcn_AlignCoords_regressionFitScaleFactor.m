function S = ...
    fcn_AlignCoords_regressionFitScaleFactor( ...
    coord_base_points, coord_xform_points, varargin)
% fcn_AlignCoords_regressionFitScaleFactor
% performs regression fitting to find the best scale factor that matches
% one set of points to another.
%
% FORMAT:
%
%  S = ...
%    fcn_AlignCoords_regressionFitScaleFactor( ...
%    coord_base_points, coord_xform_points, (fig_num))
%
% INPUTS:
%
%     coord_base_points: the points in the base coordinate system that is
%     intended to be used
%
%     coord_xform_points: the same points, in the same order, in a
%     coordinate system that needs to be changed into the base coordinate
%     system.
%
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a
%     figure to be drawn showing results.
%
% OUTPUTS:
%
%     S: the scale factor result
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_AlignCoords_regressionFitScaleFactor
% for a full test suite.
%
% This function was written on 2023_03_23 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2023_03_23 by Sean Brennan
% -- first write of function


% TO DO:
%
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plots = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 846;
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1 == flag_check_inputs
    
    % Are there the right number of inputs?
    narginchk(2,3);
    
    % Check the coord_base_points input, make sure it is '2column_of_numbers' type
    fcn_DebugTools_checkInputsToFunctions(...
        coord_base_points, '2column_of_numbers');
    
    % Check the coord_xform_points input, make sure it is
    % '2column_of_numbers' type, same number of rows as the base points
    fcn_DebugTools_checkInputsToFunctions(...
        coord_xform_points, '2column_of_numbers',length(coord_base_points(:,1)));
    
end


% Does user want to show the plots?
if 3 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plots = 1;
    end
end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%See: http://patorjk.com/software/taag/#p=display&f=Big&t=Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

% How many points are there?
Npoints = length(coord_base_points(:,1));

% Find index pairs
% These are the combinations of all point indicies that are unique
index_pairs = nchoosek(1:Npoints,2);

% Find distances
coord_base_values_squared = (coord_base_points(index_pairs(:,1),:) - coord_base_points(index_pairs(:,2),:)).^2;
coord_base_distances_squared = sum(coord_base_values_squared,2);

coord_xform_values_squared = (coord_xform_points(index_pairs(:,1),:) - coord_xform_points(index_pairs(:,2),:)).^2;
coord_xform_distances_squared = sum(coord_xform_values_squared,2);

% Perform regression to find scaling
% This is of form y = x*M
% with y being Nx1 vector, x being Nx1 vector, M is 1x1 scaling parameter
% (S^2)
y = coord_base_distances_squared.^0.5;
x = coord_xform_distances_squared.^0.5;
M = (y'*x)\(y'*y);
S = 1/M;

%§
%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if flag_do_plots

    fprintf(1,'Scale is: %.2f\n',S);

    figure(fig_num);
    clf;
    hold on
    
    % Create the inverse transform
    theta = 0;
    tx = 0;
    ty = 0;
    order_string = 's';
    T = fcn_AlignCoords_generate2DTransformMatrix( ...
        1/S, theta, tx, ty, order_string);
    
    % Before
    subplot(1,2,1);
    hold on;
    axis equal;
    grid on;
    
    % Plot the results
    Npoints = length(coord_base_points(:,1));
    
    % Connect the points with lines
    for ith_point = 1:Npoints
        plot([coord_base_points(ith_point,1) coord_xform_points(ith_point,1)],[coord_base_points(ith_point,2) coord_xform_points(ith_point,2)],'b-');
    end
    
    % Show the before/after
    plot(coord_base_points(:,1),coord_base_points(:,2),'r.-','MarkerSize',10,'LineWidth',5);
    plot(coord_xform_points(:,1),coord_xform_points(:,2),'b.-','MarkerSize',10,'LineWidth',5);
    
    title('Before');
    
    % after
    subplot(1,2,2);
    hold on;
    axis equal;
    grid on;
    
    normalized_coord_xform_points = [coord_xform_points ones(Npoints,1)];
    moved_points = (T*normalized_coord_xform_points')';

    
    % Plot the results
    Npoints = length(coord_base_points(:,1));
    
    % Connect the points with lines
    for ith_point = 1:Npoints
        plot([coord_base_points(ith_point,1) coord_xform_points(ith_point,1)],[coord_base_points(ith_point,2) coord_xform_points(ith_point,2)],'b-');
    end
    
    % Connect the points with lines
    for ith_point = 1:Npoints
        plot([moved_points(ith_point,1) coord_xform_points(ith_point,1)],[moved_points(ith_point,2) coord_xform_points(ith_point,2)],'g-');
    end
    
    % Show the before/after
    plot(coord_base_points(:,1),coord_base_points(:,2),'r.-','MarkerSize',10,'LineWidth',5);
    plot(coord_xform_points(:,1),coord_xform_points(:,2),'b.-','MarkerSize',10,'LineWidth',5);
    plot(moved_points(:,1),moved_points(:,2),'g.-','MarkerSize',10,'LineWidth',5);
    
    title('After');
    
    
    
    
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function



%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
