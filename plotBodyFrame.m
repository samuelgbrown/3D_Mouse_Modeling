function plotBodyFrame(M, scale)
% This function will plot the body frame represented by the configuration
% M.

% Perform some horrible bastardizations of mathematics to extract the unit
% vectors of the frame represented by M
shift = M(1:3, [4 4 4]);
M = (scale*M(1:3, 1:3) + shift)';
TPlot = permute(cat(3, M, shift'), [3 2 1]);

% Plot the frame
a = gca;
hold(a, 'on');
c = {'r', 'g', 'b'};
for i = 1:size(TPlot, 3)
    plot3(TPlot(:, 1, i), TPlot(:, 2, i), TPlot(:, 3, i), c{i});
end

% quiver3(M(1, [4 4 4]), M(2, [4 4 4]), M(3, [4 4 4]), M(1, 1:3), M(2, 1:3), M(3, 1:3), scale)