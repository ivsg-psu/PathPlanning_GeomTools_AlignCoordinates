function [T,err] = ...
    fcn_AlignCoords_fitAffineXform(coord_base_points, coord_xform_points, varargin)
% fcn_AlignCoords_fitAffineXform
% performs regression fitting to find the affine transform that matches
% one point set another. The affine transform is not technically a
% coordinate rotatin, but can be any transformation of points in homogenous
% form. Specifically, the 2-D corodinate transform of ortho coordinates
% from one frame to another with scale S, rotation theta, and translations
% tx and ty, given by:
%
% X' = [S*cos(theta) -S*sin(theta)  tx] * | x | 
% Y' = [S*sin(theta)  S*cos(theta)  ty]   | y |
% 1  = [      0             0       1 ]   | 1 |
%
% is a special case of the affine tranformation:
%
% X' = [M11  M12  M13] * | x | 
% Y' = [M21  M22  M23]   | y |
% 1  = [0     0    1 ]   | 1 |
%
% One can note here: fitting the affine transform requires 6 measurements
% and thus 3 points in (x,y) to fully fit, whereas the
% translation-rotation-scaling only requires 4 measurements and thus only 2
% points. One can also observe that the generic nature of the affine
% transform means that, with more degrees of freedom to fit, the affine fit will
% be better, usually, than the 2D coordinate fitting as long as the number of
% points is much larger than 3. Similarly, the affine fit may be worse if the
% number of points is only slightly larger than or equal to 3.
%
% Note also: the affine transform does NOT generate a true rotation and
% thus can add artifacts if this is used to represent true 2-D rotations of
% coordinate systems, primarily because M11 may not equal M22 and M12 may
% not equal -M21. In general, if it is known that coordinates are truly
% orthogonal in both the base and transformed representations, then the 2D
% coordinate fitting should be used. But, if the transform is not truly a
% coordinate scaling, rotation, and translation, then the affine fitting
% method of this function should be used.
%
% To solve for the affine matrix in regression form, the above form can be
% rewritten in regression form as:
%
% X  = [x   y   1   0   0   0]  * | M11 |
% Y  = [0   0   0   x   y   1]    | M12 |
%                                 | M13 |
%                                 | M21 |
%                                 | M22 |
%                                 | M23 |
%
% which is of the linear regressor form: y = X_reg*m assuming y = [X; Y]
% and X_reg is given by the matrix above, and m = [M11; M12; (etc) M23].
%
% This has the solution: m = (X_reg'*X_reg)\(X_reg'*y);
%
% FORMAT:
%
%  [T, err] = ...
%    fcn_AlignCoords_fitAffineXform( ...
%    coord_base_points, coord_xform_points, fig_num))
%
% INPUTS:
%
%     coord_base_points: the points in the base coordinate system that is
%     intended to be used, typically Nx2 or Nx3, but can be any dimensional
%     space.
%
%     coord_xform_points: the same points, in the same order, in a
%     coordinate system that needs to be changed into the base coordinate
%     system. Must be of row length N where N is the length of
%     coord_base_points.
%
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a
%     figure to be drawn showing results.
%
% OUTPUTS:
%
%     T: the transform matrix
%
%     err: the point-by-point error in the transform
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_AlignCoords_fitAffineXform
% for a full test suite.
%
% This function was written on 2023_03_23 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2023_03_23 to 2023_03_29, by Sean Brennan
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


% Set up problem:
Npoints = length(coord_base_points(:,1));
colOfOne = ones(Npoints,1);
colOfZero = zeros(Npoints,1);

y = [coord_xform_points(:,1); coord_xform_points(:,2)];
X_reg = [coord_base_points(:,1) coord_base_points(:,2) colOfOne colOfZero                 colOfZero           colOfZero;
     colOfZero               colOfZero            colOfZero coord_base_points(:,1) coord_base_points(:,2) colOfOne];
m = (X_reg'*X_reg)\(X_reg'*y);
T = [m(1:3)'; m(4:6)'; 0 0 1];

% Compare input points and output points
% start_points = coord_base_points;
% moved_points = (T_calculated*start_points')';

% Normalize the coord_xform_points
normalized_coord_xform_points = [coord_xform_points(:,1:2) colOfOne]; 

% Apply the pseudo-inverse to move points back to their source
moved_points = (T\normalized_coord_xform_points')';

% Calculate error between original base points and the points that are
% moved back, using vector sum form:
err = sum((moved_points(:,1:2)-coord_base_points(:,1:2)).^2,2).^0.5; 


% 
% fprintf(1,'Comparing points: \n');
% disp([moved_points coord_xform_points]);
% 
% figure(fig_num);
% clf
% hold on
% grid on;
% axis equal;
% 
% plot(coord_xform_points(:,1),coord_xform_points(:,2),'r.-','LineWidth',3,'MarkerSize',20);
% plot(coord_base_points(:,1),coord_base_points(:,2),'b.-','LineWidth',5,'MarkerSize',30);
% plot(fixed_points(:,1),fixed_points(:,2),'go-','LineWidth',3,'MarkerSize',15);
% 
% legend('Moved points','Target points','Moved points, Transformed');
% m = 2;
% Npoints = length(coord_base_points(:,1));
% title(['Affine transform with ' int2str(m) ' dimensions, ' int2str(Npoints) ' points, Maximum distance error: ' ...
%     num2str(max(err),6)]);
% title('Demonstration of Affine transform regression');
% 
% 
% fprintf(1,'Using fcn_AlignCoords_fit2DAffine:\n');
% fprintf(1,'Results of fitting entire transform matrix, T:\n');
% fprintf(1,'T true: \n');
% disp(T);
% fprintf(1,'T calculated: \n');
% disp(T_calculated);

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
    grid on;
    axis equal;
       
 
    
    
    figure(fig_num);
    clf;
    hold on
    
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
    plot(moved_points(:,1),moved_points(:,2),'g.-','MarkerSize',10,'LineWidth',3);
    
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
