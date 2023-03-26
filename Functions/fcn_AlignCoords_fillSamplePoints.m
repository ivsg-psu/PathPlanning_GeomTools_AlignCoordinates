function normalized_start_points = ...
    fcn_AlignCoords_fillSamplePoints(varargin)
% fcn_AlignCoords_fillSamplePoints
% fills in sample points in homogenous form, useful for testing
% transformations. Optional inputs allows users to specify the test point
% group, and/or plot the test points.
% 
% FORMAT:
% 
%  normalized_start_points = ...
%    fcn_AlignCoords_fillSamplePoints((type_of_points),(fig_num))
% 
% INPUTS:
% 
%     (none)
% 
%     (optional inputs)
%
%     type_of_points: an integer specifying the type of points to produce:
%         1: The big dipper points (default)
%         2: A random polytope with 7 points, random radii about 4, offset to (6,6)
%         3: A square of edge length 2, with bottom left at origin
%         4: Random normal
%         Note: an empty matrix, [], dafaults to 1
%
%     fig_num: any number that acts as a figure number output, causing a 
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
% See the script: script_test_fcn_AlignCoords_fillSamplePoints
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
    narginchk(0,2);
    %
    %     % Check the S input, make sure it is '1column_of_numbers' type, 1 row
    %     fcn_DebugTools_checkInputsToFunctions(...
    %         S, '1column_of_numbers',1);
    %
    %     % Check the theta input, make sure it is '1column_of_numbers' type, 1 row
    %     fcn_DebugTools_checkInputsToFunctions(...
    %         theta, '1column_of_numbers',1);
    %
    %     % Check the tx input, make sure it is '1column_of_numbers' type, 1 row
    %     fcn_DebugTools_checkInputsToFunctions(...
    %         tx, '1column_of_numbers',1);
    %
    %     % Check the ty input, make sure it is '1column_of_numbers' type, 1 row
    %     fcn_DebugTools_checkInputsToFunctions(...
    %         ty, '1column_of_numbers',1);
    %
    %     % Check the order_string input, make sure it is 'string' type
    %     fcn_DebugTools_checkInputsToFunctions(...
    %         order_string, '_of_chars');
       
end

type_of_points = 1;
% Does user want to spedify the type of points?
if 1 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        type_of_points = temp;
    end
end


% Does user want to show the plots?
if 2 == nargin
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

%         2: A random polytope
%         3: A square
%         4: Random normal

if type_of_points>4
    warning('Point types up to 4 supported. Defaulting to 1.');
    type_of_points = 1;
end

switch type_of_points
    case 1  % The big dipper - obtained from pixels in Photoshop of BigDipper_small image
        start_points = [461 357; 501 274; 555 247; 620 209; 680 231; 750 140; 698 83]; 
        
        % Image coordinates have y inverted and shifted. Need to fix
        start_points(:,2) = 540 - start_points(:,2); 
        
        title_string = 'Point type 1: the big dipper';
    case 2  % A random polytope with 7 points, random radii about 4, offset to (6,6)
        Npoints = 7;
        radii_mean = 4; 
        offset = [6 6];
        
        % Calculate a random angular spacing
        temp_point_spacing = rand(Npoints+1,1);
        temp_cum_sum = cumsum(temp_point_spacing);
        random_angular_spacing = temp_cum_sum/temp_cum_sum(end);
        random_angular_spacing = random_angular_spacing(1:end-1);
        angles = 2*pi*random_angular_spacing;
        radii = radii_mean + randn(Npoints,1);
        
        % Use sines and cosines to get xy form
        start_points = [radii.*cos(angles) radii.*sin(angles)]+ones(Npoints,1)*offset;
        
         % Close off the polygon
        start_points = [start_points;start_points(1,:)];
        
        title_string = 'Point type 2: random offset polytope';
        
    case 3  % A square of edge length 2, with bottom left at origin
        start_points = [0 0; 2 0; 2 2; 0 2; 0 0];
        title_string = 'Point type 3: square of edge 2';

    case 4  % A random normal distribution with 7 points, offset to (6,6)
        offset = [6,6];
        Npoints = 7;
        start_points = randn(Npoints,2) + ones(Npoints,1)*offset;

        title_string = 'Point type 4: random normal';

    otherwise
        error('Entered point type case that should not be accessible!');
end
        % % Calculate the square



% Add coordinate axes lines?
% start_points = [start_points;nan(1,2);[0 0; 2 0]; nan(1,2); [0 0; 0 2]];

% Convert to nomralized form
normalized_start_points = [start_points ones(length(start_points(:,1)),1)];



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
    
    figure(fig_num);
    hold on
    grid on
    axis square
    
    
    % Show the before/after
    plot(normalized_start_points(:,1),normalized_start_points(:,2),'r-','MarkerSize',10,'LineWidth',5);
    plot(normalized_start_points(:,1),normalized_start_points(:,2),'b.','MarkerSize',30,'LineWidth',5);
    title(title_string);

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
