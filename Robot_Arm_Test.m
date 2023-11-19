% What I would like to be able to do:
% Given some force f in the sensor's body frame coordinates (b-frame),
% translate the force to the arm's global coordinates (s-frame).
%
% Based on this force defined in the space frame, the kinematics (model as
% a system driven by a scaled force, dv/dt = f(v)), the location/normal of
% the nearest gate surface, and the gate force threshold, calculate the dv
% term for the next timestep, in the space frame

% Define s frame as that a hand would have when held in front of you
%   x - forward
%   y - left
%   z - up
%
% Define b frame as that a palm-up hand would have when held in front of
% you, with configuration affixed to the end of the second link
%   x - pointing from wrist to elbow
%   z - aligned along positive rotation axis of elbow
%
% Define joint 1 as the shoulder joint that rotates about s_z
% Define joint 2 as the shoulder joint that rotates about -s_y
% Define joint 3 as the elbow joint that rotates about -s_y (at 0 frame)
%
% Define the 0 frame to be both links laid flat along the -x axis
%
%   M = T_sb - [1 0 0 -2L]
%              [0 0 -1 0 ]
%              [0 1 0  0 ]
%              [0 0 0  1 ]
%
% Joint 1
%   S1 - [ 0]
%        [ 0]
%        [ 1]
%        [ 0]
%        [ 0]
%        [ 0]
%
% Joint 2
%   S2 - [ 0]
%        [-1]
%        [ 0]
%        [ 0]
%        [ 0]
%        [ 0]
%
% Joint 3
%   S3 - [ 0]
%        [-1]
%        [ 0]
%        [ 0]
%        [ 0]
%        [ L]
%
% Physical model:
% component of dv along M-gate_plane line must be less than or equal to the length of the M-gate_plane line, unless f_external > f_gate
% 	M-gate_plane line: line from getLoc(M) perpendicular to the x-y plane of M_G
% 	roughly, dv along line is: max(0, f_external - f_gate)
%
% (current trajectory places it at-or-below-gate)*min(-X, f_gate)
% 	min = returns the vector with the minimum amplitude
% 	X = M-gate_plane component of f_external
%
%   M_G_z = getUnitVec(M_G, 3)
%   M_G_loc = getLoc(M_G)
%   M_loc = getLoc(M)
%   unit(v) = v./|v|
%   distToPlaneAlongVec(loc, plane_frame, vec) = (dot(getUnitVec(plane_frame), loc) + getLoc(plane_frame))./unit(vec)
%   minVec(a, b) = returns vector with smallest magnitude
%   maxVec(a, b) = returns vector with largest magnitude
%
% Modes
% Driven mode - Arm may move on its own to keep the actuator away from the
% gate
% dv = f_scale.*f_ext_s + (dot(M_G_z, M_loc) < 0).*(M_G_z.*f_g) + (abs(dot(M_G_z, M_loc)) < tol).*minVec(-1*f_ext_s along M_G_z, M_G_z.*f_g);
%           Extra step in BOTH cases: if (dot(M_G_z, M_loc).*dot(M_G_z,
%           M_loc_new) < 0) && (dot(M_G_z, M_loc) > dot(M_G_z, M_loc_new))
%           (where M_loc_new = M_loc + dv.*dt), then must constrain
%           dv_new.*dt = dp_new = dp.*(distToPlaneAlongVec(M_loc, M_G,
%           dp)/|dp|)
%
% Resistive mode - Arm will only resist motion to keep the actuator away
% from the gate
% dv = f_scale.*f_ext_s + (dot(M_G_z, M_loc) <= 0).*minVec(-1*f_ext_s along M_G_z, M_G_z.*f_g)
%

%% Trivial case of single axis
% Assume that we have a single force acting on a 2-link arm
%
% Links are defined as unit length in the positive |y| direction
% Force is of unit magnitude, in the positive |z| direction
%       Produces a moment of magnitude 2 about the positive |x| axis

% Define any useful functions
% (https://danieltakeshi.github.io/2018/01/11/twists-and-exponential-coordinates/#:~:text=The%20matrix%20exponential%20of%20a,respect%20to%20the%20same%20frame.)
skew = @(x) [0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0];
matExpTwist = @(S, r) [expm(skew(S(1:3)*r)) ((eye(3) - expm(skew(S(1:3)*r)))*cross(S(1:3), S(4:6)) + (S(1:3)*S(1:3)'*S(4:6)*r));zeros(1,3) 1]; % Calculate the matrix exponential of the twist described by a rotation r about screw axis S
getRot = @(M) M(1:3, 1:3); % Get the rotation of the frame
getUnitVec = @(M, i) M(1:3, i)'; % Get the ith rotated (not translated) unit vector of the frame M (x: i == 1, y: i == 2, z: i == 3)
getLoc = @(M) M(1:3, 4)'; % Get the location of the frame
genConfig = @(x, y, z, rx, ry, rz) [(rotz(rz)*roty(ry)*rotx(rx)) [x;y;z];zeros(1,3) 1]; % Generate a configuration from a location and Euler angles
unit = @(v) v/norm(v);
distToPlaneAlongVec = @(loc_vec, plane_frame, traj_vec) dot(getUnitVec(plane_frame, 3), loc_vec - getLoc(plane_frame))*(1/dot(-getUnitVec(plane_frame, 3), unit(traj_vec)));
ternary = @(test, yes, no) test*yes + (~test)*no;
minVec = @(a, b) ternary(norm(a) < norm(b), a, b);
maxVec = @(a, b) ternary(norm(a) > norm(b), a, b);

% Define features of the simulation
t0 = 0; % Initial time
tf = 5; % Final time
dt = .05; % Time step (s) = 50ms (TODO: Should change this to being variable per frame)
tVec = t0:dt:tf; % Vector of timestamps
numT = length(tVec);
allM = nan(4, 4, numT); % Make space to record each M frame
allDv = nan(numT, 3); % Make space to record each dv
allDt = nan(numT, 1); % Make space to record each dt value
doPlot = true;
plotScale = 1;

% Define the physical system
L = 150; % Length of each link
M = genConfig(-2*L, 0, 0, 90, 0, 0); % Body frame, initialized (zero-configuration)
S1 = [0;0;1;0;0;0];
S2 = [0;-1;0;0;0;0];
S3 = [0;-1;0;0;0;L];

f_g = 1.5; % The threshold force of the gate
isDriven = false; % Which type force feedback are we going to use (true = driven mode, false = resistive mode)
tol = 1e-3; % Tolerance for being "on" the boundary (only used in driven mode)
f_scale = 50; % Input force scaling (lbs/(mm/s))

% Define the input force as a piecewise constant function
t_f_ext = [0 2.5 3 3.5  4  4.5 5];
y_f_ext = [3 6 0  0   0   4  0];
% f_ext_scalar = 2; % 4.17
f_ext_s = @(t) interp1(t_f_ext, y_f_ext, t, 'previous')*[.5 -.2 .3]; % Force on actuator: Forward, slightly left and more slightly up

% Things that will change when a new gate is sent! (Can be computed on the
% fly if we're SERIOUSLY memory lim...)
M_G = genConfig(0, 0, 20, -20 + 180, 0, -40); % Gate frame (gate surface is defined as x-y plane of this frame)
M_G_z = getUnitVec(M_G, 3);
M_G_loc = getLoc(M_G);

% Initialize some simulation-level parameters
tLast = tVec(1);
allM(:, :, 1) = M;
allDv(1, :) = 0;
allDt(1) = 0;

%% Start the simulation
for tInd = 2:numT
    % Perform simulation-level calculations (NOT USED IN FINAL PROGRAM)
    t = tVec(tInd);
    f_external_b = f_ext_s(t)*getRot(M); % Get the force in body-space
    
    % Set some parameters
    outsideGate = false; % To track if we cross into the gate boundary in driven mode
    
    % TODO: Based on the angles from the actuator potentiometers, calculate
    % our current rotation
    % Can use simple method from Wikipedia
    
    % TODO: Based on the forces from the actuator load cells, calculate our
    % external force f_external_b
    % Each acuator just becomes one component
    
    % Perform some prior calculations
    thisDt = t - tLast; % Calculate the timestep
    f_external_s = f_external_b*getRot(M)'; % Find the external force in the space frame by using the actuator's known frame
    M_Loc_prior = getLoc(M); % Get the location of the actuator frame
    M_Loc_prior_perp_bound = M_Loc_prior - M_G_loc; % A vector between the boundary plane and M_Loc_prior, that is perpendicular to the boundary
    boundDist_prior = dot(M_G_z, M_Loc_prior_perp_bound); % Find the distance to the boundary
    
    % Depending on which mode we're in, and where we are relative to the
    % boundary, do different things
    if isDriven
        % We are in driven mode
        
        % Check where we are relative to the boundary
        if abs(boundDist_prior) < tol
            % Note: Previously had a todo here that my external force was
            % incorrectly defined, but I think I fixed that...
            % We are "on" the gate boundary
            f_external_opposing = -dot(M_G_z, f_external_s);
            boundaryForce = min(norm(f_external_opposing), f_g); % Might technically not need the "abs", because f_external_opposing SHOULD always be positive...maybe...
        elseif boundDist_prior < 0
            % We are inside of the gated volume
            boundaryForce = f_g;
        else
            % We are outside of the gated volume
            boundaryForce = 0;
            outsideGate = true;
        end
    else
        % We are in resistive mode
        if boundDist_prior <= tol
            f_external_opposing = -dot(M_G_z, f_external_s)*M_G_z;
            boundaryForce = min(norm(f_external_opposing), f_g); % Might technically not need the "abs", because f_external_opposing SHOULD always be positive...maybe...
        else
            boundaryForce = 0;
            outsideGate = true;
        end
    end
    
    % Add the boundary force to dv
    dv = f_scale*(f_external_s + boundaryForce*M_G_z);
    
    % Calculate the new location
    M_Loc_post = M_Loc_prior + dv*thisDt;
    M_Loc_post_perp_bound = M_Loc_post - M_G_loc; % A vector between the boundary plane and M_Loc_prior, that is perpendicular to the boundary
    boundDist_post = dot(M_G_z, M_Loc_post_perp_bound); % Find the distance to the boundary
    
    % Perform the post-scale step, to correct for boundary entry/departure
    
    %     if outsideGate % Only can be true if in driven mode
    crossingBoundary = boundDist_prior*boundDist_post < -tol; % If the signs of the distance to the boundary are opposite, we have crossed the boundary
    if crossingBoundary
        % Figure out where we encountered the gate
        dp = dv*thisDt;
        scaleFactor = distToPlaneAlongVec(M_Loc_prior, M_G, dp)/norm(dp); % Location along dp where the gate was encountered (where 0 is the beginning and 1 is the end)
        
        if boundDist_post < 0
            % We will be entering the boundary this frame.  Check the size of
            % the projection of -f_external_s onto M_G_z to see what to do
            if norm(-dot(M_G_z, f_external_s)) < f_g
                % If the external force is less than the gate force, then
                % we should not let the actuator through!
                dp = dp*scaleFactor;
                dv = dp/thisDt; % TODO: Not really needed, except for visualization
            else
                % If the external force is greater than the gate force,
                % calculate the final deflection as a function of the
                % amount of time the actuator has spent in the boundary
                dv = f_scale*(f_external_s + (1 - scaleFactor)*f_g*M_G_z);
                dp = dv*thisDt;
            end
        else
            % We will be leaving the boundary this frame.  Check the size of
            % the projection of -f_external_s onto M_G_z to see what to do
            
            % Place the actuator onto the boundary
            dp = dp*scaleFactor;
            
            % Determine if the force is pointing away from the gate,
            % because it can't re-enter the gate, as it would not have left
            % if it could
            forcePointingAway = dot(f_external_s, M_G_z) > 0;
            
            % If we need to take the remaining force into account...
            if forcePointingAway
                % Add the force magnitude, scaled by the amount of time it
                % is outside of the boundary
                dp = dp + f_scale*(1 - scaleFactor)*f_external_s;
            end
            
            dv = dp/thisDt; % TODO: Not really needed, except for visualization
        end
        
        % Recalculate M_Loc_post
        M_Loc_post = M_Loc_prior + dp;
        boundDist_final = dot(M_G_z, M_Loc_post - M_G_loc); % TODO: Not really needed, except for debugging
    end
    %     end
    
    % Finally, save the values we need
    M(1:3, 4) = M_Loc_post; % Set the new location of the actuator frame
    tLast = t;
    
    % Do simulation-level calculations
    allM(:, :, tInd) = M;
    allDv(tInd, :) = dv;
    allDt(tInd) = thisDt;
end

%% Plot the results
if doPlot
    h = figure;
    a = axes;
    hold(a, 'on');
    %     axis square;
    axis equal;
    
    %     xlim([-2.5*L, .5*L]);
    %     ylim([-1.5*L, 1.5*L]);
    %     xlim([-1.5*L, 1.5*L]);
    
    for tInd = 1:(numT - 1)
        % For each time-step, plot the body frame and the dv vector
        thisM = allM(:, :, tInd);
        thisDv = allDv(tInd + 1, :); % Use the NEXT frame's dv, as that is the change that will be done on this time step's body frame
        thisDt = allDt(tInd);
        thisDp = thisDv*thisDt;
        
        plotBodyFrame(thisM, plotScale);
        quiver3(thisM(1, 4), thisM(2, 4), thisM(3, 4), thisDp(1), thisDp(2), thisDp(3), plotScale);
    end
    
    plotBodyFrame(eye(4), plotScale*10);
    plotSurfaceFrame(M_G, plotScale*10); % Plot the gate
end


%% Notes
% Rot matrices
%   Represent orientation
%       R_sb = orientation of b relative to s
%   Change reference: R_sc = R_sb*R_bc
%   Rotate vector: p' = R*p
%       Change frame of force:
%       f_sb = f_sa*R_sb
%   Rotate frame
%       Have frame c expressed in frame s, with rotation R
%           R_sc' = R*R_sc
%               Rotate c relative to s, expressed in s
%           R_sc'' = R_sc*R
%               Rotate c relative to c, expressed in s
%   R commonly indicates rotation of body frame relative to s frame
% Homogenous transformation matrix T, rot matrix R, position relative to
% space frame p
%   T = [R p;0(x3) 1]
%   Represent configuration (p is location of b relative to s)
%       T_sb = [R_sb p;0 0 0 1]
%       T_bs = [R_sb' -R_sb'p;0 0 0 1]
%   Change reference frame (same as rot matrices)
%       For point: [p_s;1] = T_sb * [p_b;1]
%   Displace point/frame
%       Left-multiply: Transformation (rot then trans) is expressed in
%       terms of first frame, produces b'
%       Right-multiply: Transformation (trans then rot) is expressed in
%       terms of second frame, produces b''
% Wrench F, moment vec m, force vec f, moment arm vec r
%   m = r x f
%   F = [m;f]
% Twist v, power p
%   v'F = p