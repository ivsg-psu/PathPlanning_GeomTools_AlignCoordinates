
%% Introduction to and Purpose of the Code
% This is the explanation of the code that can be found by running
%       script_demo_AlignCoordinates.m
% This is a script to demonstrate the functions within the AlignCoordinates
% code library. This code repo is typically located at:
%   https://github.com/ivsg-psu/PathPlanning_GeomTools_AlignCoordinates
%
% If you have questions or comments, please contact Sean Brennan at
% sbrennan@psu.edu
%
% The purpose of the code is to align two coordinate systems together.

% Useful links:
% https://web.cse.ohio-state.edu/~shen.94/681/Site/Slides_files/transformation_review.pdf



%% Revision History:
%      2023_03_27: - sbrennan@psu.edu
%      -- created a demo script of core debug utilities


%% Prep workspace
clc % Clear the console
close all % Close all figures

%% Dependencies and Setup of the Code
% The code requires several other libraries to work, namely the following
%%
% 
% * DebugTools - the repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools
% * PathClassLibrary - the repo can be found at: https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary
% 
% Each is automatically installed in a folder called "Utilities" under the
% root folder, for example:
%  ./Utilities/DebugTools/ ,
%  ./Utilities/PathClassLibrary/ .
% 
% For ease of transfer, zip files of the directories used - without the
% .git repo information, to keep them small - are referenced and are NOT
% included in this repo. These dependencies are to open code repos, and
% this code accesses these and downloads, if needed, the appropriate
% releases.



% USE THE FOLLOWING CODE TO ALLOW MANUAL INSTALLS OF LIBRARIES (IF PRIVATE)
% 
% The following code checks to see if the folders flag has been
% initialized, and if not, it calls the DebugTools function that loads the
% path variables. It then loads the PathClassLibrary functions as well.
% Note that the PathClass Library also has sub-utilities that are included.

% if ~exist('flag_AlignCoords_Folders_Initialized','var')
%     
%     % add necessary directories for function creation utility 
%     %(special case because folders not added yet)
%     debug_utility_folder = fullfile(pwd, 'Utilities', 'DebugTools');
%     debug_utility_function_folder = fullfile(pwd, 'Utilities', 'DebugTools','Functions');
%     debug_utility_folder_inclusion_script = fullfile(pwd, 'Utilities', 'DebugTools','Functions','fcn_DebugTools_addSubdirectoriesToPath.m');
%     if(exist(debug_utility_folder_inclusion_script,'file'))
%         current_location = pwd;
%         cd(debug_utility_function_folder);
%         fcn_DebugTools_addSubdirectoriesToPath(debug_utility_folder,{'Functions','Data'});
%         cd(current_location);
%     else % Throw an error?
%         error('The necessary utilities are not found. Please add them (see README.md) and run again.');
%     end
%     
%     % Now can add the Path Class Library automatically
%     utility_folder_PathClassLibrary = fullfile(pwd, 'Utilities', 'PathClassLibrary');
%     fcn_DebugTools_addSubdirectoriesToPath(utility_folder_PathClassLibrary,{'Functions','Utilities'});
%     
%     % utility_folder_GetUserInputPath = fullfile(pwd, 'Utilities', 'GetUserInputPath');
%     % fcn_DebugTools_addSubdirectoriesToPath(utility_folder_GetUserInputPath,{'Functions','Utilities'});
% 
%     % Now can add all the other utilities automatically
%     folder_LapsClassLibrary = fullfile(pwd);
%     fcn_DebugTools_addSubdirectoriesToPath(folder_LapsClassLibrary,{'Functions'});
% 
%     % set a flag so we do not have to do this again
%     flag_AlignCoords_Folders_Initialized = 1;
% end

% Use the following code to install public libraries

% List what libraries we need, and where to find the codes for each
clear library_name library_folders library_url

ith_library = 1;
library_name{ith_library}    = 'DebugTools_v2023_01_29';
library_folders{ith_library} = {'Functions','Data'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_29.zip?raw=true';

% ith_library = ith_library+1;
% library_name{ith_library}    = 'PathClass_v2023_02_01';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary/blob/main/Releases/PathClass_v2023_02_01.zip?raw=true';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'GetUserInputPath_v2023_02_01';
% library_folders{ith_library} = {''};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_GetUserInputPath/blob/main/Releases/GetUserInputPath_v2023_02_01.zip?raw=true';

% Do we need to set up the work space?
if ~exist('flag_AlignCoords_Folders_Initialized','var')

    % Reset all flags for installs to empty
    clear global FLAG*

    fprintf(1,'Installing utilities necessary for code ...\n');

    % Dependencies and Setup of the Code
    % This code depends on several other libraries of codes that contain
    % commonly used functions. We check to see if these libraries are installed
    % into our "Utilities" folder, and if not, we install them and then set a
    % flag to not install them again.
    
    % Set up libraries
    for ith_library = 1:length(library_name)
        dependency_name = library_name{ith_library};
        dependency_subfolders = library_folders{ith_library};
        dependency_url = library_url{ith_library};

        fprintf(1,'\tAdding library: %s ...',dependency_name);
        fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url);
        clear dependency_name dependency_subfolders dependency_url
        fprintf(1,'Done.\n');
    end

    % Set dependencies for this project specifically
    fcn_DebugTools_addSubdirectoriesToPath(pwd,{'Functions','Data'});
    
    disp('Done setting up libraries, adding each to MATLAB path, and adding current repo folders to path.');
    
end

%% Create a set of test situations
% 1) Exact match
% 2) Noisy but exact match
% 3) Very noisy match

% We choose homogenous coordinates of the form below which includes:
% translation via (tx, ty), rotation by angle (q), and uniform scaling of
% the axes (S):
%
% X' = [S*cos(q)  -sin(q)   tx ] * | x | 
% Y' = [sin(q)   S*cos(q)   ty ]   | y |
% 1  = [0             0      1  ]   | 1 |

% Fill in a sample transform matrix
S = 2;
rotation_angle = 70*pi/180;
tx = 2;
ty = 7;
Tscale = [...
         S                   0                 0; 
         0                   S                 0; 
         0                   0                 1 
    ];
Trotate = [...
     cos(rotation_angle)  -sin(rotation_angle) 0; 
      sin(rotation_angle)  cos(rotation_angle) 0; 
         0                   0                 1 
    ];

Ttranslate = [...
         1                   0                 tx; 
         0                   1                 ty; 
         0                   0                 1 
    ];

Tsrt = Tscale*Trotate*Ttranslate;
Trts = Trotate*Ttranslate*Tscale;
T = Trts;
%T = Tsrt;

% Fill in some sample points, in homogenous form

% % Calculate the square
% start_points = [0 0; 2 0; 2 2; 0 2];

% Calculate the weird polygon
Npoints = 7;
temp_point_spacing = rand(Npoints+1,1);
temp_cum_sum = cumsum(temp_point_spacing);
normalized_temp_cum_sum = temp_cum_sum/temp_cum_sum(end);
normalized_temp_cum_sum = normalized_temp_cum_sum(1:end-1);
angles = 2*pi*normalized_temp_cum_sum;
radii = 4+rand(Npoints,1);
start_points = [radii.*cos(angles) radii.*sin(angles)]+ones(Npoints,1)*[6 6];

% Close off the polygon
start_points = [start_points;start_points(1,:)];

% Add coordinate axes lines?
% start_points = [start_points;nan(1,2);[0 0; 2 0]; nan(1,2); [0 0; 0 2]];

% Convert to nomralized form
normalized_start_points = [start_points ones(length(start_points(:,1)),1)];

% Apply the transform
moved_points = (T*normalized_start_points')';

% Plot the results
fig_num = 2384;
figure(fig_num);
clf;
hold on;
grid on;
axis equal

Npoints = length(start_points(:,1));
for ith_point = 1:Npoints
    plot([start_points(ith_point,1) moved_points(ith_point,1)],[start_points(ith_point,2) moved_points(ith_point,2)],'-');
end

plot(start_points(:,1),start_points(:,2),'r.-','MarkerSize',10,'LineWidth',5);
plot(moved_points(:,1),moved_points(:,2),'b.-','MarkerSize',10,'LineWidth',5);

%% Set up solution for finding the scale factor
coord_base_points = start_points;
coord_xform_points = moved_points;

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
fprintf(1,'Scale is: %.2f\n',S);

%% Set up solution for finding the rotation

% Find start and end points
starting_edges = coord_base_points(index_pairs(:,1));
ending_edges = coord_base_points(index_pairs(:,2));

URHERE

%% Perform regression to find scaling
% This is of form y = x*M
% with y being Nx1 vector, x being Nx1 vector, M is 1x1 scaling parameter
% (S^2)
y = coord_base_distances_squared.^0.5;
x = coord_xform_distances_squared.^0.5;
M = (y'*x)\(y'*y);
S = 1/M;
fprintf(1,'Scale is: %.2f\n',S);



%% Generic affine transform 
% This has 4 parameters, and thus can be fit with minimum 2 points



% Note that the translation of ortho coordinates from one frame to another
% is a special case of the affine tranformation
% X' = [M11  M12  M13] * | x | 
% Y' = [M21  M22  M23]   | y |
% 1  = [0     0    1 ]   | 1 |
% But note that the rotation-translation-scaling transform requires 3 points to fully fit,
% whereas the translation-rotation-scaling only requires 2!






function fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url, varargin)
%% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES - MATLAB package installer from URL
%
% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES installs code packages that are
% specified by a URL pointing to a zip file into a default local subfolder,
% "Utilities", under the root folder. It also adds either the package
% subfoder or any specified sub-subfolders to the MATLAB path.
%
% If the Utilities folder does not exist, it is created.
% 
% If the specified code package folder and all subfolders already exist,
% the package is not installed. Otherwise, the folders are created as
% needed, and the package is installed.
% 
% If one does not wish to put these codes in different directories, the
% function can be easily modified with strings specifying the
% desired install location.
% 
% For path creation, if the "DebugTools" package is being installed, the
% code installs the package, then shifts temporarily into the package to
% complete the path definitions for MATLAB. If the DebugTools is not
% already installed, an error is thrown as these tools are needed for the
% path creation.
% 
% Finally, the code sets a global flag to indicate that the folders are
% initialized so that, in this session, if the code is called again the
% folders will not be installed. This global flag can be overwritten by an
% optional flag input.
%
% FORMAT:
%
%      fcn_DebugTools_installDependencies(...
%           dependency_name, ...
%           dependency_subfolders, ...
%           dependency_url)
%
% INPUTS:
%
%      dependency_name: the name given to the subfolder in the Utilities
%      directory for the package install
%
%      dependency_subfolders: in addition to the package subfoder, a list
%      of any specified sub-subfolders to the MATLAB path. Leave blank to
%      add only the package subfolder to the path. See the example below.
%
%      dependency_url: the URL pointing to the code package.
%
%      (OPTIONAL INPUTS)
%      flag_force_creation: if any value other than zero, forces the
%      install to occur even if the global flag is set.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      This code will automatically get dependent files from the internet,
%      but of course this requires an internet connection. If the
%      DebugTools are being installed, it does not require any other
%      functions. But for other packages, it uses the following from the
%      DebugTools library: fcn_DebugTools_addSubdirectoriesToPath
%
% EXAMPLES:
%
% % Define the name of subfolder to be created in "Utilities" subfolder
% dependency_name = 'DebugTools_v2023_01_18';
%
% % Define sub-subfolders that are in the code package that also need to be
% % added to the MATLAB path after install; the package install subfolder
% % is NOT added to path. OR: Leave empty ({}) to only add 
% % the subfolder path without any sub-subfolder path additions. 
% dependency_subfolders = {'Functions','Data'};
%
% % Define a universal resource locator (URL) pointing to the zip file to
% % install. For example, here is the zip file location to the Debugtools
% % package on GitHub:
% dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_18.zip?raw=true';
%
% % Call the function to do the install
% fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
%
% This function was written on 2023_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_23:
% -- wrote the code originally

% TO DO
% -- Add input argument checking

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
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

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(3,4);
end

%% Set the global variable - need this for input checking
% Create a variable name for our flag. Stylistically, global variables are
% usually all caps.
flag_varname = upper(cat(2,'flag_',dependency_name,'_Folders_Initialized'));

% Make the variable global
eval(sprintf('global %s',flag_varname));

if nargin==4
    if varargin{1}
        eval(sprintf('clear global %s',flag_varname));
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist(flag_varname,'var') || isempty(eval(flag_varname))
    % Save the root directory, so we can get back to it after some of the
    % operations below. We use the Print Working Directory command (pwd) to
    % do this. Note: this command is from Unix/Linux world, but is so
    % useful that MATLAB made their own!
    root_directory_name = pwd;

    % Does the directory "Utilities" exist?
    utilities_folder_name = fullfile(root_directory_name,'Utilities');
    if ~exist(utilities_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(root_directory_name,'Utilities');

        % Did it work?
        if ~success_flag
            error('Unable to make the Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The Utilities directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
        end

    end

    % Does the directory for the dependency folder exist?
    dependency_folder_name = fullfile(root_directory_name,'Utilities',dependency_name);
    if ~exist(dependency_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(utilities_folder_name,dependency_name);

        % Did it work?
        if ~success_flag
            error('Unable to make the dependency directory: %s. Reason: %s with message ID: %s\n',dependency_name, error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The %s directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',dependency_name, error_message, message_ID);
        end

    end

    % Do the subfolders exist?
    flag_allFoldersThere = 1;
    if isempty(dependency_subfolders)
        flag_allFoldersThere = 0;
    else
        for ith_folder = 1:length(dependency_subfolders)
            subfolder_name = dependency_subfolders{ith_folder};
            
            % Create the entire path
            subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
            
            % Check if the folder and file exists that is typically created when
            % unzipping.
            if ~exist(subfunction_folder,'dir')
                flag_allFoldersThere = 0;
            end
        end
    end

    % Do we need to unzip the files?
    if flag_allFoldersThere==0
        % Files do not exist yet - try unzipping them.
        save_file_name = tempname(root_directory_name);
        zip_file_name = websave(save_file_name,dependency_url);
        % CANT GET THIS TO WORK --> unzip(zip_file_url, debugTools_folder_name);

        % Is the file there?
        if ~exist(zip_file_name,'file')
            error('The zip file: %s for dependency: %s did not download correctly. This is usually because permissions are restricted on the current directory. Check the code install (see README.md) and try again.\n',zip_file_name, dependency_name);
        end

        % Try unzipping
        unzip(zip_file_name, dependency_folder_name);

        % Did this work?
        flag_allFoldersThere = 1;
        if ~isempty(dependency_subfolders)
            for ith_folder = 1:length(dependency_subfolders)
                subfolder_name = dependency_subfolders{ith_folder};
                
                % Create the entire path
                subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
                
                % Check if the folder and file exists that is typically created when
                % unzipping.
                if ~exist(subfunction_folder,'dir')
                    flag_allFoldersThere = 0;
                end
            end
        end
        
        if flag_allFoldersThere==0
            error('The necessary dependency: %s has an error in install, or error performing an unzip operation. Check the code install (see README.md) and try again.\n',dependency_name);
        else
            % Clean up the zip file
            delete(zip_file_name);
        end

    end


    % For path creation, if the "DebugTools" package is being installed, the
    % code installs the package, then shifts temporarily into the package to
    % complete the path definitions for MATLAB. If the DebugTools is not
    % already installed, an error is thrown as these tools are needed for the
    % path creation.
    %
    % In other words: DebugTools is a special case because folders not
    % added yet, and we use DebugTools for adding the other directories
    if strcmp(dependency_name(1:10),'DebugTools')
        debugTools_function_folder = fullfile(root_directory_name, 'Utilities', dependency_name,'Functions');

        % Move into the folder, run the function, and move back
        cd(debugTools_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        cd(root_directory_name);
    else
        try
            fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        catch
            error('Package installer requires DebugTools package to be installed first. Please install that before installing this package');
        end
    end


    % Finally, the code sets a global flag to indicate that the folders are
    % initialized.  Check this using a command "exist", which takes a
    % character string (the name inside the '' marks, and a type string -
    % in this case 'var') and checks if a variable ('var') exists in matlab
    % that has the same name as the string. The ~ in front of exist says to
    % do the opposite. So the following command basically means: if the
    % variable named 'flag_CodeX_Folders_Initialized' does NOT exist in the
    % workspace, run the code in the if statement. If we look at the bottom
    % of the if statement, we fill in that variable. That way, the next
    % time the code is run - assuming the if statement ran to the end -
    % this section of code will NOT be run twice.

    eval(sprintf('%s = 1;',flag_varname));
end

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

    % Nothing to do!



end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function fcn_DebugTools_installDependencies