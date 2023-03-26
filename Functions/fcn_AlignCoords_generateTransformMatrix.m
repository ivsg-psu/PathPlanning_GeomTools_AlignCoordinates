function T = ...
    fcn_AlignCoords_generateTransformMatrix( ...
    S, theta, tx, ty, order_string...
    )
% fcn_AlignCoords_generateTransformMatrix
% fills in a transform matrix given scaling, angle, and translation values.
% 
% FORMAT:
% 
%    [ ...
%    T
%    ] = ...
%    fcn_AlignCoords_generateTransformMatrix( ...
%    seed_points, ...
%    V, ...
%    C, ...
%    AABB, ...
%    stretch, ...
%    (fig_num) ...
%    )
% 
% INPUTS:
% 
%     seed_points: the list of seed points in [x y] format, where x and y 
%     are columns
% 
%     V: the V matrix resulting from Voronoi calculations
% 
%     C: the C matrix resulting from Voronoi calculations
%
%     AABB: the axis-aligned bounding box, in format of 
%     [xmin ymin xmax ymax], wherein the resulting polytopes must be
%     bounded.
% 
%     stretch: the factor to stretch the polytopes after they are
%     calculated, in format of [xmagnification ymagnification]
% 
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a 
%     figure to be drawn showing results.

Tsrt = Tscale*Trotate*Ttranslate
% 
% 
% OUTPUTS:
% 
%     polytopes: the resulting polytopes after converting to polytope form.
% 
% 
% DEPENDENCIES:
% 
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_cropPolytopeToRange
%     fcn_MapGen_fillPolytopeFieldsFromVertices
% 
% 
% EXAMPLES:
% 
% See the script: script_test_fcn_MapGen_generatePolysFromVoronoiAABBWithTiling
% for a full test suite.
% 
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

% 
% REVISION HISTORY:
% 
% 2021_07_02 by Sean Brennan
% -- first write of function
% 2021_07_30 by Sean Brennan
% -- fixed errors due to corners being omitted
 
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
    narginchk(3,4);

    % Check the seed_points input, make sure it is '2column_of_numbers' type
    fcn_DebugTools_checkInputsToFunctions(...
        seed_points, '2column_of_numbers');
 
    % Check the AABB input, make sure it is '4column_of_numbers' type
    fcn_DebugTools_checkInputsToFunctions(...
        AABB, '4column_of_numbers',1);
    
    % Check the stretch input, make sure it is '2column_of_numbers' type
    fcn_DebugTools_checkInputsToFunctions(...
        stretch, '2column_of_numbers',1);
       
end


% Does user want to show the plots?
if 4 == nargin
    fig_num = varargin{end};
    if ~isempty(fig_num)
        figure(fig_num);
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
    clf;
    hold on
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);
    
    % plot the polytopes
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',2);
    
    % plot all vertices
    plot(all_vertices(:,2),all_vertices(:,3),'c','Linewidth',1);
    
    
    % plot the seed points in red
    plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);
    
    % plot the means in black
    temp = zeros(length(polytopes),2);
    for ith_poly = 1:length(polytopes)
        temp(ith_poly,:) = polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3);
    
    % number the polytopes at seed points
    for ith_poly = 1:length(polytopes)
        text_location = seed_points(ith_poly,:);
        text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
    end
        
    %     % number the polytopes at means
    %     for ith_poly = 1:length(polytopes)
    %         text_location = polytopes(ith_poly).mean;
    %         text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
    %     end
    
    % plot the connections between the polytope neighbors
    if 1==0
        % Clean up and sort the vertices so that we can associate neighbors
        all_vertices_no_nan = all_vertices(~isnan(all_vertices(:,1)),:);
        sorted_all_vertices = sortrows(all_vertices_no_nan,[2 3]);
        
        % Remove repeats
        sorted_all_vertices = unique(sorted_all_vertices,'rows','stable');
        
        % Remove infinities
        sorted_all_vertices = sorted_all_vertices(~isinf(sorted_all_vertices(:,2)));
        
        Nrealvertices = floor(length(sorted_all_vertices(:,1))/3);
        data = zeros(Nrealvertices*6,2);
        for ith_poly = 1:Nrealvertices
            row_offset = (ith_poly-1)*3;
            neighbors = sorted_all_vertices(row_offset+1:row_offset+3,1);
            
            for jth_neighbor = 2:length(neighbors)
                neigh_offset = (ith_poly-1)*6 + ((jth_neighbor-2)*3);
                data(neigh_offset+1:neigh_offset+3,:) = [seed_points(neighbors(1),:); seed_points(neighbors(jth_neighbor),:); nan nan];
            end
        end
        plot(data(:,1),data(:,2),'-','Linewidth',0.5);
    end
    
    
    %% Show a detailed step-by-step process behind construction of obstacle map 
    if 1==0
        % using fcn_MapGen_generatePolysFromVoronoiAABB
        fig_num = 1010;
        figure(fig_num);
        clf;
        
        % Calculate the scale
        scale = max(AABB,[],'all') - min(AABB,[],'all');
        new_axis = [AABB(1)-1.5*scale AABB(3)+1.5*scale AABB(2)-1.5*scale AABB(4)+1.5*scale];


        %% plot the seed points in red
        subplot(2,3,1);
        
        plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);
        
        
        % number the polytopes at seed points
        for ith_poly = 1:length(polytopes)
            text_location = seed_points(ith_poly,:);
            text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
        end
        axis(new_axis);
        title('Seed points');
        
        %% plot the tiled seed points
        subplot(2,3,2);
        cla;
        plot(tiled_original_seed_points(:,1),tiled_original_seed_points(:,2),'b.','Markersize',20);
        hold on;
        plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);

        axis(new_axis);
        title('Tiled Seed points');
        
        
        %% PLOT THE VORONOI lines with the points
        subplot(2,3,3);
        
        % Start the loop to calculate all the vertices
        num_poly = size(tiled_original_seed_points,1);
        voronoi_polytopes(num_poly) = ...
            struct(...
            'seed_point',[],...
            'vertices',[],...
            'AABB',[],...
            'xv',[],...
            'yv',[],...
            'distances',[],...
            'mean',[],...
            'area',[],...
            'max_radius',[],...
            'min_radius',[],...
            'mean_radius',[],...
            'radii',[],...
            'cost',[],...
            'parent_poly_id',[]);
        
        Npolys = length(voronoi_polytopes);
        Nvertices_per_poly = 20; % Maximum estimate
        Nvertices_per_map = Npolys*Nvertices_per_poly;
        all_voronoi_vertices = nan(Nvertices_per_map,3);

        for ith_poly = 1:length(tiled_original_seed_points(:,1))
             
            % Fill in seed_point
            voronoi_polytopes(ith_poly).seed_point = tiled_original_seed_points(ith_poly);
            
            % Get the verticies
            vertices_open = V(C{ith_poly},:);
            
            % Append results to close off the vector loop
            vertices = [vertices_open; vertices_open(1,:)]; % Close off the vertices
            
            % Save verticies to this polytope
            voronoi_polytopes(ith_poly).vertices = vertices;
                        
            % Save verticies to the all_verticies array
            Nvertices = length(vertices(:,1));
            if Nvertices>Nvertices_per_poly
                error('Need to resize the number of allowable vertices');
            else
                row_offset = (ith_poly-1)*Nvertices_per_poly;
                all_voronoi_vertices(row_offset+1:row_offset+Nvertices,1) = ith_poly;
                all_voronoi_vertices(row_offset+1:row_offset+Nvertices,2:3) = vertices;
            end
        end

        
        % plot the tiled seed points in blue
        plot(tiled_original_seed_points(:,1),tiled_original_seed_points(:,2),'b.','Markersize',20);
        hold on;
        
        % plot all vertices
        plot(all_voronoi_vertices(:,2),all_voronoi_vertices(:,3),'c','Linewidth',1);
        
        axis(new_axis);
        title('Voronoi boundaries');
        
        
        %% PLOT THE VORONOI lines for just the middle
        subplot(2,3,4);
        
        % plot the polytopes on current axis
        fcn_MapGen_plotPolytopes(polytopes,gca,'b',2);
        hold on;
        
        
        % plot the seed points in red
        plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);
        
        axis(new_axis);
        title('Extract middle, tilable polytopes')
        
        
        %% show that it's tilable
        subplot(2,3,5);
                
        % plot the polytopes
        fcn_MapGen_plotPolytopes(polytopes,gca,'b',5);
        hold on;
        
        xscale = AABB(3)-AABB(1);
        yscale = AABB(4)-AABB(2);
        
        shifted_left_up = fcn_INTERNAL_shift_polys(polytopes,-xscale,yscale);
        fcn_MapGen_plotPolytopes(shifted_left_up,gca,'-',2);
              
        shifted_left = fcn_INTERNAL_shift_polys(polytopes,-xscale,0);
        fcn_MapGen_plotPolytopes(shifted_left,gca,'-',2);
                
        shifted_left_down = fcn_INTERNAL_shift_polys(polytopes,-xscale,-yscale);
        fcn_MapGen_plotPolytopes(shifted_left_down,gca,'-',2);
                
        shifted_up = fcn_INTERNAL_shift_polys(polytopes,0,yscale);
        fcn_MapGen_plotPolytopes(shifted_up,gca,'-',2);
                
        shifted_down = fcn_INTERNAL_shift_polys(polytopes,0,-yscale);
        fcn_MapGen_plotPolytopes(shifted_down,gca,'-',2);
        
        shifted_right_up = fcn_INTERNAL_shift_polys(polytopes,xscale,yscale);
        fcn_MapGen_plotPolytopes(shifted_right_up,gca,'-',2);
        
        shifted_right = fcn_INTERNAL_shift_polys(polytopes,xscale,0);
        fcn_MapGen_plotPolytopes(shifted_right,gca,'-',2);
        
        shifted_right_down = fcn_INTERNAL_shift_polys(polytopes,xscale,-yscale);
        fcn_MapGen_plotPolytopes(shifted_right_down,gca,'-',2);
        
        axis(new_axis);
        title('Polytopes can tile');
        
        %% plot the shrink to edge
        subplot(2,3,6);
        
        des_gap_size = 0.05;
        
        shrunk_polytopes=...
            fcn_MapGen_polytopesShrinkFromEdges(...
            polytopes,des_gap_size);
        
        % plot the shrunk polytopes
        fcn_MapGen_plotPolytopes(shrunk_polytopes,gca,'r',2);
        hold on;
        
        
        shifted_left_up = fcn_INTERNAL_shift_polys(shrunk_polytopes,-xscale,yscale);
        fcn_MapGen_plotPolytopes(shifted_left_up,gca,'b-',2);
              
        shifted_left = fcn_INTERNAL_shift_polys(shrunk_polytopes,-xscale,0);
        fcn_MapGen_plotPolytopes(shifted_left,gca,'b-',2);
                
        shifted_left_down = fcn_INTERNAL_shift_polys(shrunk_polytopes,-xscale,-yscale);
        fcn_MapGen_plotPolytopes(shifted_left_down,gca,'b-',2);
                
        shifted_up = fcn_INTERNAL_shift_polys(shrunk_polytopes,0,yscale);
        fcn_MapGen_plotPolytopes(shifted_up,gca,'b-',2);
                
        shifted_down = fcn_INTERNAL_shift_polys(shrunk_polytopes,0,-yscale);
        fcn_MapGen_plotPolytopes(shifted_down,gca,'b-',2);
        
        shifted_right_up = fcn_INTERNAL_shift_polys(shrunk_polytopes,xscale,yscale);
        fcn_MapGen_plotPolytopes(shifted_right_up,gca,'b-',2);
        
        shifted_right = fcn_INTERNAL_shift_polys(shrunk_polytopes,xscale,0);
        fcn_MapGen_plotPolytopes(shifted_right,gca,'b-',2);
        
        shifted_right_down = fcn_INTERNAL_shift_polys(shrunk_polytopes,xscale,-yscale);
        fcn_MapGen_plotPolytopes(shifted_right_down,gca,'b-',2);
               
        axis(new_axis);
        title('Shrunk from edge polytopes, with tiling');
    end
    
    if 1==0
        %% Confirm the tiling
        % if it is a true tiling, then the parts that stick out, e.g. values
        % greater than the AABB boundary, should match points that are in the
        % interior - in other words, there is a seamless connection between the two
        % areas.
        
        % Find the parts that stick out
        all_points = all_vertices(~isnan(all_vertices(:,1)),2:3);
        all_points_shifted = all_points - AABB(1:2);
        rounded_points = mod(all_points_shifted,1);
        x_indices_stick_out = find(all_points_shifted(:,1)~=rounded_points(:,1));
        y_indices_stick_out = find(all_points_shifted(:,2)~=rounded_points(:,2));
        xy_indices_stick_out = intersect(x_indices_stick_out,y_indices_stick_out);
        
        % Find the points that are entirely within the AABB
        x_indices_inside = find(all_points_shifted(:,1)==rounded_points(:,1));
        y_indices_inside = find(all_points_shifted(:,2)==rounded_points(:,2));
        xy_indices_inside = intersect(x_indices_inside,y_indices_inside);
        points_inside = all_points(xy_indices_inside,:);
        
        % DEBUG AREA
        tiling_fig_num = 12312;
        figure(tiling_fig_num);
        clf;
        hold on
        
        scale = max(AABB,[],'all') - min(AABB,[],'all');
        new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
        axis(new_axis);
        
        % plot the polytopes
        fcn_MapGen_plotPolytopes(polytopes,tiling_fig_num,'b',2);
        
        % plot all vertices
        plot(all_points(:,1),all_points(:,2),'c.','Linewidth',1);
        
        % Plot all the points that stick out
        plot(all_points(x_indices_stick_out,1),all_points(x_indices_stick_out,2),'r.','Markersize',20);
        plot(all_points(y_indices_stick_out,1),all_points(y_indices_stick_out,2),'m.','Markersize',15);
        plot(all_points(xy_indices_stick_out,1),all_points(xy_indices_stick_out,2),'k.','Markersize',10);
        legend('Polytopes','All points','X-data out of bounds','Y-data out of bounds','Both X and Y out of bounds');
        
        % For each of the "stick out" points, check to see that there's an inside
        % point that corresponds to the same location, but inside.
        
        legend off;
        
        % Start with x
        error_inside_to_outside_x = nan(length(x_indices_stick_out),1);
        for ith_outside_point = 1:length(x_indices_stick_out)
            current_outside_point = all_points(x_indices_stick_out(ith_outside_point),:);
            rounded_current_outside_point = mod(current_outside_point,1);
            
            % Use a distance metric to find closest actual point
            tolerance = 1E-10;
            distances_to_outside_point = sum((rounded_current_outside_point-points_inside).^2,2).^0.5;
            [min_distance,~] = min(distances_to_outside_point);
            index_min = distances_to_outside_point<tolerance;
            closest_points = points_inside(index_min,:);
            
            % Show the results
            plot(current_outside_point(:,1),current_outside_point(:,2),'ro','Markersize',20);
            plot(closest_points(:,1),closest_points(:,2),'ro','Markersize',20);
            error_inside_to_outside_x(ith_outside_point)= min_distance;
        end
        
        % Now with y
        error_inside_to_outside_y = nan(length(y_indices_stick_out),1);
        for ith_outside_point = 1:length(y_indices_stick_out)
            current_outside_point = all_points(y_indices_stick_out(ith_outside_point),:);
            rounded_current_outside_point = mod(current_outside_point,1);
            
            % Use a distance metric to find closest actual point
            tolerance = 1E-10;
            distances_to_outside_point = sum((rounded_current_outside_point-points_inside).^2,2).^0.5;
            [min_distance,~] = min(distances_to_outside_point);
            index_min = distances_to_outside_point<tolerance;
            closest_points = points_inside(index_min,:);
            
            % Show the results
            plot(current_outside_point(:,1),current_outside_point(:,2),'go','Markersize',20);
            plot(closest_points(:,1),closest_points(:,2),'go','Markersize',20);
            error_inside_to_outside_y(ith_outside_point)= min_distance;
        end
        figure(1234);
        clf;
        hold on;
        plot(error_inside_to_outside_x);
        plot(error_inside_to_outside_y);
        legend('X errors', 'Y errors');
    end
    
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
