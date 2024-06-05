function out = sphereplot(data, referencePoints, HRTFtarget)
figure;
% Unique target positions
targetPositions = unique(data(:, 1));

% Prepare the tiled layout
tiledlayout('flow');

% Define a colormap
colors = lines(2); % Two colors for the two stimuli

% Store plot handles for common legend
plotHandles = [];

% Loop through each target position
for i = 1:numel(targetPositions)
    % Filter data for the current target position
    currentPosition = targetPositions(i);
    filteredData = data(data(:, 1) == currentPosition, :);
    out{i} = filteredData;

    % Next tile
    nexttile;
    hold on;

    % Plot the hemisphere guide with circular gradient color
    [X, Y, Z] = sphere(50);
    azimuth = atan2(Y, X); % Azimuth angles
    elevation = asin(Z);   % Elevation angles
    distance = sqrt(azimuth.^2 + elevation.^2);
    distance = distance / max(distance(:));
    distance = 1 - distance;
    surf(X, Y, Z, distance, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

    % Plot the data for both stimuli
    for j = 0:1
        stimulusData = filteredData(filteredData(:, 2) == j, :);
        azimuthRad = abs(deg2rad(stimulusData(:, 3)));
        elevationRad = deg2rad(stimulusData(:, 4));

        % Calculate x, y, z coordinates on the unit sphere
        x = cos(elevationRad) .* cos(azimuthRad);
        y = cos(elevationRad) .* sin(azimuthRad);
        z = sin(elevationRad);

        % Plot the data points and store the plot handle
        h = scatter3(x, y, z, 50, colors(j+1, :), 'filled', 'DisplayName', ['Stimulus ' num2str(j)]);
        plotHandles = [plotHandles, h]; % Collect plot handles for legend
    end

    % Plot the reference point for the current target position
    refPoint = referencePoints(referencePoints(:, 1) == currentPosition, 2:3);
    refAzimuthRad = deg2rad(refPoint(1));
    refElevationRad = deg2rad(refPoint(2));
    refX = cos(refElevationRad) * cos(refAzimuthRad);
    refY = cos(refElevationRad) * sin(refAzimuthRad);
    refZ = sin(refElevationRad);

    h = scatter3(refX, refY, refZ, 100, 'k', 'filled', 'DisplayName', 'Reference Point');
    plotHandles = [plotHandles, h]; % Collect plot handles for legend
    text(refX, refY, refZ, '  Ref', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

    % Set plot labels and title
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(sprintf('Ref Position (Az: %d°, El: %d°)', refPoint(1), refPoint(2)));
    grid on;
    axis equal;
    view(3); % Ensure 3D view

    hold off;
end

% Set uniform axis limits
allAxes = findall(gcf, 'Type', 'axes');
for k = 1:numel(allAxes)
    set(allAxes(k), 'XLim', [-1 1], 'YLim', [-1 1], 'ZLim', [-1 1]);
end

% Overall plot title
sgtitle(sprintf('[%s] Unit Sphere Plot of Azimuth and Elevation by Target Position and Stimulus', HRTFtarget), 'Interpreter', 'none');

% Create a common legend
legend(plotHandles(1:3), 'Stimulus 0', 'Stimulus 1', 'Reference Point', 'Location', 'best');
end