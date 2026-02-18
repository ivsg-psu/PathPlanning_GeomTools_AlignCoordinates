function T = ...
    fcn_AlignCoords_generate2DTransformMatrix( ...
    S, theta, tx, ty, order_string, varargin)
% fcn_AlignCoords_generate2DTransformMatrix
% fills in a transform matrix given scaling, angle, and translation values.
% It uses the homogenous form of the transformation in 2D, allowing for
% scaling, rotation, and translation given by:
%
% Tscale = [...
%          S                   0                 0; 
%          0                   S                 0; 
%          0                   0                 1 
%     ];
% 
% 
% Trotate = [...
%      cos(theta)  -sin(theta) 0; 
%       sin(theta)  cos(theta) 0; 
%          0         0         1 
%     ];
% 
% Ttranslate = [...
%          1                   0                 tx; 
%          0                   1                 ty; 
%          0                   0                 1 
%     ];
%
% The transformation of homogenous from one coordinate system to another is
% achieved by: moved_point = (T*start_spoints')' if the start points are
% Nx3, or by moved_point = (T*start_spoints) if the start points are 3xN.
% 
% The inverse transformation can be found by the pseudo-inverse: 
% inverse_points = (T\normalized_coord_xform_points')';
%
% For example, assume that coord_base_points is filled with data with N
% rows and 2 columns. The following code demonstrates creation of a
% transform, moving points, then using the pseudo-inverse tranformation to
% move them back:
%
%         figure; % Set up a new figure
%         hold on; % Hold the plots
%         
%         % Define a transformation
%         S = 2;
%         theta = 45*pi/180;
%         tx = 3;
%         ty = 7;
%         order_string = 'tsr';
%         T = fcn_AlignCoords_generate2DTransformMatrix(S, theta, tx, ty, order_string);
%         
%         Npoints = 10; % How many points should we demo with?
%         coord_base_points = rand(Npoints,2); % Set up some test points
%         plot(coord_base_points(:,1),coord_base_points(:,2),'b.','Markersize',10); % Plot the result
%         
%         % Convert to homogenous form in prep for transform
%         colOfOne = ones(Npoints,1); % Define column of 1's
%         normalized_coord_base_points = [coord_base_points(:,1:2) colOfOne];  % Make homogenous
%         normalized_coord_xform_points = (T*normalized_coord_base_points')';
%         
%         % Plot the results
%         plot(normalized_coord_xform_points(:,1),normalized_coord_xform_points(:,2),'r.','Markersize',10); % Plot the result
%         
%         % Invert them back
%         moved_back_points = (T\normalized_coord_xform_points')';
%         
%         % Plot the results, and show that the moved back points, given by green
%         % circles, land back at the original points, shown in blue
%         plot(moved_back_points(:,1),moved_back_points(:,2),'go','Markersize',10); % Plot the result
%
% 
% FORMAT:
% 
%  T = ...
%     fcn_AlignCoords_generateTransformMatrix( ...
%     S, theta, tx, ty, order_string,
%     (figNum) ...
%     )
% 
% INPUTS:
% 
%     S: the scaling factor
% 
%     theta: the angle of rotation
% 
%     tx: the translation in x
% 
%     tx: the translation in y
%
%     order_string: a string of the form 'srt' that indicates the order of
%     operations to create matrix T. For example, 
%     'srt' creates a matrix: Tsrt = Tscale*Trotate*Ttranslate
%     'trs' creates a matrix: Ttrs = Ttranslate*Trotate*Tscale
%     't' creates a matrix: Ttrs = Ttranslate
%     (etc)
% 
%     (optional inputs)
%
%     figNum: any number that acts as a figure number output, causing a 
%     figure to be drawn showing results.
% 
% OUTPUTS:
% 
%     T: a transformation matrix in homogenous coordinates
% 
% 
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
% 
% EXAMPLES:
% 
% See the script: script_test_fcn_AlignCoords_generate2DTransformMatrix
% for a full test suite.
% 
% This function was written on 2023_03_23 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

% 
% REVISION HISTORY:
% 
% 2023_03_23 by Sean Brennan
% - first write of function

 
% TO DO:
% 
% - fill in to-do items here.

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
    narginchk(5,6);

    % Check the S input, make sure it is '1column_of_numbers' type, 1 row
    fcn_DebugTools_checkInputsToFunctions(...
        S, '1column_of_numbers',1);
 
    % Check the theta input, make sure it is '1column_of_numbers' type, 1 row
    fcn_DebugTools_checkInputsToFunctions(...
        theta, '1column_of_numbers',1);
    
    % Check the tx input, make sure it is '1column_of_numbers' type, 1 row
    fcn_DebugTools_checkInputsToFunctions(...
        tx, '1column_of_numbers',1);
    
    % Check the ty input, make sure it is '1column_of_numbers' type, 1 row
    fcn_DebugTools_checkInputsToFunctions(...
        ty, '1column_of_numbers',1);
    
    % Check the order_string input, make sure it is 'string' type
    fcn_DebugTools_checkInputsToFunctions(...
        order_string, '_of_chars');
       
end


% Does user want to show the plots?
if 6 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        figNum = temp;
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
Tscale = [...
         S                   0                 0; 
         0                   S                 0; 
         0                   0                 1 
    ];


Trotate = [...
     cos(theta)  -sin(theta) 0; 
      sin(theta)  cos(theta) 0; 
         0         0         1 
    ];

Ttranslate = [...
         1                   0                 tx; 
         0                   1                 ty; 
         0                   0                 1 
    ];
T = eye(3);

for ith_character = 1:length(order_string)
    current_character = order_string(ith_character);
    switch lower(current_character)
        case 'r'
            T = T*Trotate;
        case {'s'}
            T = T*Tscale;
        case 't'
            T = T*Ttranslate;
        otherwise
            error('Unknown rotation type: %s', current_character);
    end
end


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
    
    figure(figNum);
    clf;
    hold on
    axis equal
    grid on;
    
    % Calculate the weird polygon
    rng(23);
    start_points = fcn_AlignCoords_fillSamplePoints;
    
    % Apply the transform
    moved_points = (T*start_points')';
    
    % Plot the results
    Npoints = length(start_points(:,1));
    
    % Connect the points with lines
    for ith_point = 1:Npoints
        plot([start_points(ith_point,1) moved_points(ith_point,1)],[start_points(ith_point,2) moved_points(ith_point,2)],'-');
    end
    
    % Show the before/after
    plot(start_points(:,1),start_points(:,2),'r.-','MarkerSize',10,'LineWidth',3);
    plot(start_points(:,1),start_points(:,2),'k.','MarkerSize',20,'LineWidth',3);
    plot(moved_points(:,1),moved_points(:,2),'b.-','MarkerSize',10,'LineWidth',3);
    plot(moved_points(:,1),moved_points(:,2),'k.','MarkerSize',20,'LineWidth',3);

    title(sprintf('Results of operation: %s',upper(order_string)));
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
