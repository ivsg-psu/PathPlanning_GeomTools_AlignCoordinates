%% Introduction to and Purpose of the Code
% This is the explanation of the code that can be found by running
%       
%    script_demo_AlignCoordinates.m
% 
% This is a script to demonstrate the functions within the AlignCoordinates
% code library. This code repo is typically located at:
% 
%   https://github.com/ivsg-psu/PathPlanning_GeomTools_AlignCoordinates
%
% If you have questions or comments, please contact Sean Brennan at
% sbrennan@psu.edu
%
% The purpose of the code is to align two coordinate systems together.

% Useful links:
% https://web.cse.ohio-state.edu/~shen.94/681/Site/Slides_files/transformation_review.pdf



% REVISION HISTORY:
% 
% 2023_03_27 by Sean Brennan, sbrennan@psu.edu
% - Created a demo script of core debug utilities
%
% 2026_02_18 by Sean Brennan, sbrennan@psu.edu
% - In script_demo_AlignCoordinates
%   % * Updated headers to current standard
%   % * Added auto-installer
%
% (new release)

% TO-DO:
% - 2026_02_18 by Sean Brennan, sbrennan@psu.edu
%   * Need to update all functions and scripts to standard form


%% Make sure we are running out of root directory
st = dbstack; 
thisFile = which(st(1).file);
[filepath,name,ext] = fileparts(thisFile);
cd(filepath);

%%% START OF STANDARD INSTALLER CODE %%%%%%%%%

%% Clear paths and folders, if needed
if 1==1
    clear flag_AlignCoordinates_Folders_Initialized
end

if 1==0
    fcn_INTERNAL_clearUtilitiesFromPathAndFolders;
end

if 1==0
    % Resets all paths to factory default
    restoredefaultpath;
end

%% Install dependencies
% Define a universal resource locator (URL) pointing to the repos of
% dependencies to install. Note that DebugTools is always installed
% automatically, first, even if not listed:
clear dependencyURLs dependencySubfolders
ith_repo = 0;

% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary';
% dependencySubfolders{ith_repo} = {'Functions','Data'};
 
ith_repo = ith_repo+1;
dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_PathTools_GetUserInputPath';
dependencySubfolders{ith_repo} = {''};

% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_GeomTools_AlignCoordinates';
% dependencySubfolders{ith_repo} = {'Functions'};
% 
% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/FieldDataCollection_GPSRelatedCodes_GPSClass';
% dependencySubfolders{ith_repo} = {'Functions'};



% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/FieldDataCollection_VisualizingFieldData_PlotRoad';
% dependencySubfolders{ith_repo} = {'Functions','Data'};


% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary';
% dependencySubfolders{ith_repo} = {'Functions','testFixtures','GridMapGen'};



%% Do we need to set up the work space?
if ~exist('flag_AlignCoordinates_Folders_Initialized','var')
    
    % Clear prior global variable flags
    clear global FLAG_*

    % Navigate to the Installer directory
    currentFolder = pwd;
    cd('Installer');
    % Create a function handle
    func_handle = @fcn_DebugTools_autoInstallRepos;

    % Return to the original directory
    cd(currentFolder);

    % Call the function to do the install
    func_handle(dependencyURLs, dependencySubfolders, (0), (-1));

    % Add this function's folders to the path
    this_project_folders = {...
        'Functions','Data'};
    fcn_DebugTools_addSubdirectoriesToPath(pwd,this_project_folders)

    flag_AlignCoordinates_Folders_Initialized = 1;
end

%%% END OF STANDARD INSTALLER CODE %%%%%%%%%

%% Set environment flags for input checking in library
% These are values to set if we want to check inputs or do debugging
setenv('MATLABFLAG_ALIGNCOORDINATES_FLAG_CHECK_INPUTS','1');
setenv('MATLABFLAG_ALIGNCOORDINATES_FLAG_DO_DEBUG','0');

%% Set environment flags that define the ENU origin
% This sets the "center" of the ENU coordinate system for all plotting
% functions
% Location for Test Track base station
setenv('MATLABFLAG_PLOTROAD_REFERENCE_LATITUDE','40.86368573');
setenv('MATLABFLAG_PLOTROAD_REFERENCE_LONGITUDE','-77.83592832');
setenv('MATLABFLAG_PLOTROAD_REFERENCE_ALTITUDE','344.189');


%% Set environment flags for plotting
% These are values to set if we are forcing image alignment via Lat and Lon
% shifting, when doing geoplot. This is added because the geoplot images
% are very, very slightly off at the test track, which is confusing when
% plotting data
setenv('MATLABFLAG_PLOTROAD_ALIGNMATLABLLAPLOTTINGIMAGES_LAT','-0.0000008');
setenv('MATLABFLAG_PLOTROAD_ALIGNMATLABLLAPLOTTINGIMAGES_LON','0.0000054');

%% Check repo compliance?
if 1==0
	% Set input arguments
	repoShortName = '_AlignCoordinates_';

	% Call the function
	fcn_DebugTools_testRepoForRelease(repoShortName, (2222));
end

%% Start of Demo Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _             _            __   _____                          _____          _
%  / ____| |           | |          / _| |  __ \                        / ____|        | |
% | (___ | |_ __ _ _ __| |_    ___ | |_  | |  | | ___ _ __ ___   ___   | |     ___   __| | ___
%  \___ \| __/ _` | '__| __|  / _ \|  _| | |  | |/ _ \ '_ ` _ \ / _ \  | |    / _ \ / _` |/ _ \
%  ____) | || (_| | |  | |_  | (_) | |   | |__| |  __/ | | | | | (_) | | |___| (_) | (_| |  __/
% |_____/ \__\__,_|_|   \__|  \___/|_|   |_____/ \___|_| |_| |_|\___/   \_____\___/ \__,_|\___|
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Start%20of%20Demo%20Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Welcome to the demo code for the AlignCoordinates library!')

% We choose homogenous coordinates of the form below which includes:
% translation via (tx, ty), rotation by angle (q), and uniform scaling of
% the axes (S):
%
% X' = [S*cos(q)  -S*sin(q)   tx ] * | x | 
% Y' = [S*sin(q)   S*cos(q)   ty ]   | y |
% 1  = [0             0       1  ]   | 1 |



%% Demonstrate fcn_AlignCoords_generate2DTransformMatrix
% Fills in transformation matricies, and can demonstrate the transform on
% sample sets of points

fig_num = 1;

% Fill in a sample transform matrix
S = 2;
tx = 2;
ty = 7;
theta = -40*pi/180;

order_string = 'RTS';
T = fcn_AlignCoords_generate2DTransformMatrix( ...
    S, theta, tx, ty, order_string, fig_num);
title('Demonstration of fcn_AlignCoords_fillSamplePoints', 'Interpreter', 'none');

%% Demonstrate fcn_AlignCoords_fillSamplePoints
% Fills in some sample points, in homogenous form
fig_num = 2;
figure(fig_num);
clf;

for type_of_points = 1:4
    subplot(2,2,type_of_points);
    start_points = fcn_AlignCoords_fillSamplePoints(type_of_points, fig_num); %#ok<NASGU>
end
sgtitle('Demonstration of fcn_AlignCoords_generate2DTransformMatrix', 'Interpreter', 'none','FontSize',12);


%% Ccreate test points used for all functions that follow
start_points = fcn_AlignCoords_fillSamplePoints;
moved_points = (T*start_points')';

noise_level = 0.1;
noise_to_points = [noise_level*randn(length(start_points(:,1)),2) ones(length(start_points(:,1)),1)];

coord_base_points = start_points;
coord_xform_points = moved_points + noise_to_points;

%% Show how to apply the transform
moved_points = (T*start_points')';

% Plot the results
fig_num = 3;
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

title('Demonstration of applying transforms to homogenous points');

%% Demonstrate fcn_AlignCoords_regressionFitScaleFactor
% Shows function for finding the scale factor
fig_num = 4;
S_calculated = fcn_AlignCoords_regressionFitScaleFactor(coord_base_points(:,1:2), coord_xform_points(:,1:2), fig_num);
sgtitle('Demonstration of fcn_AlignCoords_regressionFitScaleFactor', 'Interpreter', 'none','FontSize',12);

fprintf(1,'Using fcn_AlignCoords_regressionFitScaleFactor:\n');
fprintf(1,'Results of fitting scale factor, S:\n');
fprintf(1,'S true: \n');
disp(S);

fprintf(1,'S calculated: \n');
disp(S_calculated);

%% Demonstrate fcn_AlignCoords_fitRotationKabsch
% Show solution for finding the rotation and translation

fig_num = 5;

% NOTE: this function can work in hyperdimensions or even 1D, so we have to
% carefully specify to only use the 2 columns
[R_calculated,t_rotated,err] = fcn_AlignCoords_fitRotationKabsch(coord_base_points(:,1:2), coord_xform_points(:,1:2),fig_num); % Find optimal transform
sgtitle('Demonstration of fcn_AlignCoords_fitRotationKabsch', 'Interpreter', 'none','FontSize',12);

% Calculate unrotated version of t, and in column form. Note: it's in the
% negative direction because the fitting method is from moved to base.
t_calculated = -t_rotated/R_calculated;

fprintf(1,'Using fcn_AlignCoords_fitRotationKabsch:\n');
fprintf(1,'Results of fitting rotation matrix, R:\n');
fprintf(1,'R true: \n');
R = T(1:2,1:2); % Extract rotation matrix out of the T matrix
disp(R);

fprintf(1,'R calculated: \n');
disp(R_calculated);

fprintf(1,'Results of fitting translation vector, t:\n');
fprintf(1,'t true: \n');
t = T(1:2,3); % Extract the t vector out of the T matrix
disp(t);
fprintf(1,'t calculated: \n');
disp(t_calculated);

fprintf(1,'Resulting errors:\n');
disp(err);


%
% fig_num = 6;
% figure(fig_num);
% clf;
% 
% hold on;

% P3 = bsxfun(@plus,coord_xform_points*R_calculated,t_rotated);
% plot(coord_xform_points(:,1),coord_xform_points(:,2),'r.-','LineWidth',3,'MarkerSize',20);
% plot(coord_base_points(:,1),coord_base_points(:,2),'b.-','LineWidth',3,'MarkerSize',20);
% plot(P3(:,1),P3(:,2),'ko-','LineWidth',3,'MarkerSize',10);
% 
% axis equal;
% grid on;
% legend('Moved points','Target points','Moved points Transformed');
% m = 2;
% Npoints = length(coord_base_points(:,1));
% title([int2str(m) ' dimensions, ' int2str(Npoints) ' points, Maximum RMSE: ' ...
%     num2str(max(err),6)]);


%% Demonstrate fcn_AlignCoords_fit2DCoordinates
% Solution for generic 2D coordinate matching (scaling, rotation, and
% translation)

fig_num = 7;
[T_calculated,R_calculated,S_calculated,t_calculated,err] = fcn_AlignCoords_fit2DCoordinates(coord_base_points(:,1:2), coord_xform_points(:,1:2), fig_num); % Find optimal transform
sgtitle('Demonstration of fcn_AlignCoords_fit2DCoordinates', 'Interpreter', 'none','FontSize',12);


fprintf(1,'Using fcn_AlignCoords_fit2DCoordinates:\n');
fprintf(1,'Results of fitting entire transform matrix, T:\n');
fprintf(1,'T true: \n');
disp(T);
fprintf(1,'T calculated: \n');
disp(T_calculated);

if 1 == 0 % The following assertions should work if the noise above is removed.
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
end

maxerr_2D = max(err,[],'all');


%% Generic affine transform 
% Show that the affine transform works, but it's not a pure rotation

fig_num = 8;
[T_calculated,err] = fcn_AlignCoords_fitAffineXform(coord_base_points(:,1:2), coord_xform_points(:,1:2), fig_num); % Find optimal transform
sgtitle('Demonstration of fcn_AlignCoords_fitAffineXform', 'Interpreter', 'none','FontSize',12);


fprintf(1,'Using fcn_AlignCoords_fit2DCoordinates:\n');
fprintf(1,'Results of fitting entire transform matrix, T:\n');
fprintf(1,'T true: \n');
disp(T);
fprintf(1,'T calculated: \n');
disp(T_calculated);

if 1 == 0 % The following assertions should work if the noise above is removed.
    % Is the error small?
    assert(max(err,[],'all')<1E-10);    
end

maxerr_affine = max(err,[],'all');



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