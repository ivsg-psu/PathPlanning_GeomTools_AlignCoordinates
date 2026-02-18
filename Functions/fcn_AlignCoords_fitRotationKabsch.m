function [R,t,err] = ...
    fcn_AlignCoords_fitRotationKabsch(coord_base_points, coord_xform_points, varargin)
% fcn_AlignCoords_fitRotationKabsch
% performs regression fitting to find the rotation that matches
% one set of points to another using the Kabsch algorithm.
%
% FORMAT:
%
%  [R, t, err] = ...
%    fcn_AlignCoords_fitRotationKabsch( ...
%    coord_base_points, coord_xform_points, figNum))
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
%     figNum: any number that acts as a figure number output, causing a
%     figure to be drawn showing results.
%
% OUTPUTS:
%
%     R: the rotation matrix
%
%     t: the translation from one coordinate system to another
%
%     err: the Least Root Mean Square Error (LRMS) in the transform
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_AlignCoords_fitRotationKabsch
% for a full test suite.
%
% This function was written on 2023_03_23 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu
%
% NOTE: this is just a wrapped version of: rigidtform
% https://www.mathworks.com/matlabcentral/fileexchange/61697-rigidtform
% Andrew Horchler (2023). rigidtform (https://github.com/horchler/rigidtform), GitHub. Retrieved March 25, 2023.
% *********
%
% Andrew D. Horchler, *horchler @ gmail . com*,
% [biorobots.case.edu](http://biorobots.case.edu/) Created: 12-8-16,
% Revision: 1.0, 2-17-17
%
% This version tested with Matlab 9.1.0.441655 (R2016b) Mac OS X 10.12.3
% (Build: 16D32), Java 1.7.0_75-b13 Compatibility maintained back through
% Matlab 8.3 (R2014a) &nbsp;
%
% *********
%
% Copyright &copy; 2016&ndash;2017, Andrew D. Horchler All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%  * Redistributions of source code must retain the above copyright notice,
%  this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%  notice, this list of conditions and the following disclaimer in the
%  documentation and/or other materials provided with the distribution.
%  * Neither the name of Case Western Reserve University nor the names of
%  its contributors may be used to endorse or promote products derived from
%  this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.
%
% See also:
% Ehud Schreiber (2023). Kabsch algorithm (https://www.mathworks.com/matlabcentral/fileexchange/25746-kabsch-algorithm), MATLAB Central File Exchange. Retrieved March 23, 2023.
% A copy of this appears at bottom:


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
    narginchk(2,3);
    
    % Check the coord_base_points input, make sure it is '2column_of_numbers' type
    fcn_DebugTools_checkInputsToFunctions(...
        coord_base_points, '1or9column_of_numbers');
    
    % Check the coord_xform_points input, make sure it is
    % '2column_of_numbers' type, same number of rows as the base points
    fcn_DebugTools_checkInputsToFunctions(...
        coord_xform_points, '1or9column_of_numbers',length(coord_base_points(:,1)));
    
end


% Does user want to show the plots?
if 3 == nargin
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
P2 = coord_base_points;
P1 = coord_xform_points;

[R,t,err] = rigidtform(P1,P2);

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
    
    
    moved_points = bsxfun(@plus,P1*R,t);              % Apply optimal transform to points
    maxerr = max(err);                      % Maximum RMSE
    
    % Plot points
    m = length(P1(1,:)); % What is the dimension of the points?
    n = length(P1(:,1)); % How many points?
    
    
    
    % Only do the debug plot in 2 dimensions
    if m == 2
%         plot(P1(:,1),P1(:,2),'ro',...
%             P2(:,1),P2(:,2),'bo',...
%             P3(:,1),P3(:,2),'k.');
%         
%         axis equal;
%         grid on;
%         legend('Moved points','Target points','Moved points Transformed');

        
        
        figure(figNum);
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
        
        %         title([int2str(m) ' dimensions, ' int2str(n) ' points, 
        
        title('After');
        
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

function [R,t,err] = rigidtform(P1,P2,dim)
%RIGIDTFORM  Optimal rotation and translation between corresponding points
%
%   [R,T] = RIGIDTFORM(P1,P2) returns the rotation matrix, R, and
%   translation vector, T, for the optimal rigid transform between two sets
%   of corresponding N-dimensional points, P1 and P2. Each row of the
%   equal-sized, real-valued floating point matrices P1 and P2 represents
%   the N-dimensional coordinates of a point. The outputs, R and T, are
%   N-by-N and 1-by-N arrays, respectively. The transform can be applied as
%   P1(i,:)*R+T to rotate and translate P1 to P2 or as (P2(i,:)-T)*R' to go
%   from P2 to P1.
%
%   [...] = RIGIDTFORM(P1,P2,DIM) optionally specifies which dimension the
%   points are stored in of the two datasets, P1 and P2. The default is 1,
%   i.e., each row represents the coordinates of a point. If DIM is 2, each
%   column represents the coordinates of a point and the outputs R and T
%   are transposed such that the transform can be applied as R*P1(:,i)+T to
%   rotate and translate P1 to P2.
%
%   [R,T,ERR] = RIGIDTFORM(...) returns the root mean squared error between
%   each corresponding point in the datasets P1 and P2 as a column vector.
%   If DIM is 2, ERR is transposed.
%
%   Example:
%       m = 3; n = 50; R_orig = orth(rand(m)); t_orig = rand(1,m);
%       P1 = rand(n,m); P2 = P1*R_orig+repmat(t_orig,n,1);
%       [R,t,err] = rigidtform(P1,P2); P3 = P1*R+repmat(t,n,1);
%       figure;
%       plot3(P1(:,1),P1(:,2),P1(:,3),'ro',...
%             P2(:,1),P2(:,2),P2(:,3),'bo',...
%             P3(:,1),P3(:,2),P3(:,3),'k.');
%       axis equal; grid on; title(['Maximum error: ' num2str(max(err))]);
%       legend('P1','P2','P1 Transformed');
%
%   See also: SVD, IMREGTFORM, AFFINE2D, AFFINE3D

%   An N-dimensional implementation of the "Kabsch algorithm":
%
%   [1] Wolfgang Kabsch (1976) "A solution for the best rotation to relate
%   two sets of vectors", Acta Crystallographica 32:922.
%   http://dx.doi.org/10.1107/S0567739476001873
%
%   [2] Wolfgang Kabsch (1978) "A discussion of the solution for the best
%   rotation to relate two sets of vectors", Acta Crystallographica 34:827.
%   http://dx.doi.org/doi:10.1107/S0567739478001680

%   Andrew D. Horchler, horchler @ gmail . com, Created 12-8-16
%   Revision: 1.0, 2-20-17
%   https://github.com/horchler/rigidtform/blob/master/rigidtform.m


% Check input points
if ~ismatrix(P1) || ~isfloat(P1) || ~isreal(P1) || ~all(isfinite(P1(:)))
    error('rigidtform:InvalidDataset1',...
        'First dataset must be a finite real floating point matrix.');
end
if ~ismatrix(P2) || ~isfloat(P2) || ~isreal(P2) || ~all(isfinite(P2(:)))
    error('rigidtform:InvalidDataset2',...
        'Second dataset must be a finite real floating point matrix.');
end

% Check that there are sufficient and equal number of input points
s = size(P1);
if ~isequal(s,size(P2))
    error('rigidtform:SizeMismatch',...
        'Size of input datasets must be equal.');
end

% Check dimension, transpose inputs if needed
if nargin < 3
    dim = 1;
else
    if ~isscalar(dim)
        error('rigidtform:NonscalarDimension',...
            'Dimension must be a scalar.');
    end
    if dim ~= 1 && dim ~= 2
        error('rigidtform:InvalidDimension','Dimension must be 1 or 2.');
    end
    if dim == 2
        P1 = P1';
        P2 = P2';
        s = s([2 1]);
    end
end

if s(2) > s(1) && s(2) > 2^15
    warning('rigidtform:LargeDimension',...
        ['Number of dimensions (%d) is very large. Datasets may be '...
        'transposed.'], s(2));
end

% Calculate centroids
c1 = mean(P1,1);
c2 = mean(P2,1);

% Center data, calculate H matrix, and return unitary matrices
[U,~,V] = svd(bsxfun(@minus,P1,c1)'*bsxfun(@minus,P2,c2));

% Calculate optimal rotation matrix
R = U*V';

% Correct R if needed to ensure right-handed coordinate system
if s(2) == 3 && det(R) < 0
    V(:,end) = -V(:,end);
    R = U*V';
end

% Calculate optimal translation vector
t = c2-c1*R;

% Calculate root mean squared distance error for each point, if requested
if nargout > 2
    d = P2-(P1*R+t);
    err = sqrt(sum(d.^2,2));
    
    % Check for and handle overflow
    idx = isinf(err);
    if any(idx)
        f = d(idx,1);
        for i = 2:s(2)
            f = hypot(f,d(idx,i));
        end
        err(idx) = f;
    end
end

% Return transposed outputs if needed
if dim == 2
    R = R';
    t = t';
    if nargout > 2
        err = err';
    end
end
end  % Ends rigidtform

%% kabsch algorithm
% https://www.mathworks.com/matlabcentral/fileexchange/25746-kabsch-algorithm

function[U, r, lrms] = Kabsch(P, Q, m) %#ok<DEFNU>
% Find the Least Root Mean Square distance
% between two sets of N points in D dimensions
% and the rigid transformation (i.e. translation and rotation)
% to employ in order to bring one set that close to the other,
% Using the Kabsch (1976) algorithm.
% Note that the points are paired, i.e. we know which point in one set
% should be compared to a given point in the other set.
%
% References:
% 1) Kabsch W. A solution for the best rotation to relate two sets of vectors. Acta Cryst A 1976;32:9223.
% 2) Kabsch W. A discussion of the solution for the best rotation to relate two sets of vectors. Acta Cryst A 1978;34:8278.
% 3) http://cnx.org/content/m11608/latest/
% 4) http://en.wikipedia.org/wiki/Kabsch_algorithm
%
% We slightly generalize, allowing weights given to the points.
% Those weights are determined a priori and do not depend on the distances.
%
% We work in the convention that points are column vectors;
% some use the convention where they are row vectors instead.
%
% Input  variables:
%  P : a D*N matrix where P(a,i) is the a-th coordinate of the i-th point
%      in the 1st representation
%  Q : a D*N matrix where Q(a,i) is the a-th coordinate of the i-th point
%      in the 2nd representation
%  m : (Optional) a row vector of length N giving the weights, i.e. m(i) is
%      the weight to be assigned to the deviation of the i-th point.
%      If not supplied, we take by default the unweighted (or equal weighted)
%      m(i) = 1/N.
%      The weights do not have to be normalized;
%      we divide by the sum to ensure sum_{i=1}^N m(i) = 1.
%      The weights must be non-negative with at least one positive entry.
% Output variables:
%  U : a proper orthogonal D*D matrix, representing the rotation
%  r : a D-dimensional column vector, representing the translation
%  lrms: the Least Root Mean Square
%
% Details:
% If p_i, q_i are the i-th point (as a D-dimensional column vector)
% in the two representations, i.e. p_i = P(:,i) etc., and for
% p_i' = U p_i + r      (' does not stand for transpose!)
% we have p_i' ~ q_i, that is,
% lrms = sqrt(sum_{i=1}^N m(i) (p_i' - q_i)^2)
% is the minimal rms when going over the possible U and r.
% (assuming the weights are already normalized).
%

sz1 = size(P) ;
sz2 = size(Q) ;
if (length(sz1) ~= 2 || length(sz2) ~= 2)
    error 'P and Q must be matrices' ;
end
if (any(sz1 ~= sz2))
    error 'P and Q must be of same size' ;
end
D = sz1(1) ;         % dimension of space
N = sz1(2) ;         % number of points
if (nargin >= 3)
    if (~isvector(m) || any(size(m) ~= [1 N]))
        error 'm must be a row vector of length N' ;
    end
    if (any(m < 0))
        error 'm must have non-negative entries' ;
    end
    msum = sum(m) ;
    if (msum == 0)
        error 'm must contain some positive entry' ;
    end
    m = m / msum ;     % normalize so that weights sum to 1
else                 % m not supplied - use default
    m = ones(1,N)/N ;
end
p0 = P*m' ;          % the centroid of P
q0 = Q*m' ;          % the centroid of Q
v1 = ones(1,N) ;     % row vector of N ones
P = P - p0*v1 ;      % translating P to center the origin
Q = Q - q0*v1 ;      % translating Q to center the origin
% C is a covariance matrix of the coordinates
% C = P*diag(m)*Q'
% but this is inefficient, involving an N*N matrix, while typically D << N.
% so we use another way to compute Pdm = P*diag(m),
% which is equivalent to, but more efficient than,
% Pdm = zeros(D,N) ;
% for i=1:N
% 	Pdm(:,i) = m(i)*P(:,i) ;
% end
Pdm = bsxfun(@times,m,P) ;
C = Pdm*Q' ;
%	C = P*Q' / N ;       % (for the non-weighted case)
[V,~,W] = svd(C) ;   % singular value decomposition
I = eye(D) ;
if (det(V*W') < 0)   % more numerically stable than using (det(C) < 0)
    I(D,D) = -1 ;
end
U = W*I*V' ;
r = q0 - U*p0 ;
Diff = U*P - Q ;     % P, Q already centered
%	lrms = sqrt(sum(sum(Diff.*Diff))/N) ; % (for the non-weighted case)
% to compute the lrms, we employ an efficient method, equivalent to:
% lrms = 0 ;
% for i=1:N
% 	lrms = lrms + m(i)*Diff(:,i)'*Diff(:,i) ;
% end
% lrms = sqrt(lrms) ;
lrms = sqrt(sum(sum(bsxfun(@times,m,Diff).*Diff))) ;
end % Ends Kabsch
