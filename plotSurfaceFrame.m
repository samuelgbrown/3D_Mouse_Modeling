function plotSurfaceFrame(M, normalScale)
% This function will plot the surface represented by the x-y plane of the
% configuration M.

% Set some parameters
surfScale = 1e5; % Arbitrarily large scale (so that the patch goes beyond the current edges of the axis)
xLims = xlim; % Get the limits of the current axis (so we can reset them later)
yLims = ylim;
zLims = zlim;

% zVec = M(1:3, 3);

% planeEq = @(x, y) -(x*zVec(1) + y*zVec(2) + dot(zVec, M(1:3, 4)))/zVec(3);

% xCol = surfScale*[-1 -1 1 1]';
% yCol = surfScale*[1 -1 -1 1]';
% zCol = planeEq(xCol, yCol);

% Perform some horrible bastardizations of mathematics to extract the unit
% vectors of the frame represented by M
M_scaleCrop = surfScale*M(1:3, 1:2);
M_scaleCropT = M_scaleCrop';
M_rotSurf = [M_scaleCropT;-M_scaleCropT];
shift = M(1:3, [4 4 4 4])';
allPatches = M_rotSurf + shift;

normalShift = shift(1, :);
% normalVec = [(normalScale*M(1:3, 3)' + normalShift);normalShift];
normalVec = M(1:3, 3);

% Plot the patch
hold(gca, 'on');
fill3(allPatches(:, 1), allPatches(:, 2), allPatches(:, 3), 'r', 'facealpha', .2);
quiver3(normalShift(1), normalShift(2), normalShift(3), normalVec(1), normalVec(2), normalVec(3), normalScale);
% fill3(xCol, yCol, zCol, 'r', 'faceAlpha', .2);
% plot3(normalVec(:, 1), normalVec(:, 2), normalVec(:, 3), 'r');

% Reset the axis
xlim(xLims);
ylim(yLims);
zlim(zLims);