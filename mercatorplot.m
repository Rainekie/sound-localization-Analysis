function out = mercatorplot(data, referencePoints, HRTFtarget)
figure;

% Unique target positions
targetPositions = unique(data(:, 1));

% Prepare the tiled layout
t = tiledlayout('flow');
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Define a colormap
colors = lines(2); % Two colors for the two stimuli

% Store plot handles for common legend
plotHandles = [];

% Loop through each target position
for i = 1:numel(targetPositions)
    % for i = 5
    % Filter data for the current target position
    currentPosition = targetPositions(i);
    filteredData = data(data(:, 1) == currentPosition, :);
    out(:,:,i) = filteredData;

    % Next tile
    nexttile;
    hold on;

    % Plot the data for both stimuli
    for j = 0:1
        stimulusData = filteredData(filteredData(:, 2) == j, :);
        azimuthDeg = abs(stimulusData(:, 3)); % Azimuth in degrees
        elevationDeg = stimulusData(:, 4); % Elevation in degrees
        azimuthRad = deg2rad(azimuthDeg);
        elevationRad = deg2rad(elevationDeg);

        % Apply Mercator projection
        [x, y] = mercatorProjection(azimuthRad, elevationRad);

        % Plot the data points and store the plot handle
        h = scatter(x, y, 50, colors(j+1, :), 'filled', 'DisplayName', ['Stimulus ' num2str(j)]);
        if isempty(plotHandles) || ~any(arrayfun(@(handle) isequal(handle, h), plotHandles))
            plotHandles = [plotHandles, h]; % Collect plot handles for legend
        end
    end

    % Plot the reference point for the current target position
    refPoint = referencePoints(referencePoints(:, 1) == currentPosition, 2:3);
    refAzimuthDeg = refPoint(1); % Reference azimuth in degrees
    refElevationDeg = refPoint(2); % Reference elevation in degrees
    refAzimuthRad = deg2rad(refAzimuthDeg);
    refElevationRad = deg2rad(refElevationDeg);

    % Apply Mercator projection
    [refX, refY] = mercatorProjection(refAzimuthRad, refElevationRad);

    h = scatter(refX, refY, 100, 'k', 'filled', 'DisplayName', 'Reference Point');
    if isempty(plotHandles) || ~any(arrayfun(@(handle) isequal(handle, h), plotHandles))
        plotHandles = [plotHandles, h]; % Collect plot handles for legend
    end
    text(refX, refY, '  Ref', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

    % Set plot labels and title
    xlabel('Azimuth (degrees)');
    ylabel('Elevation (degrees)');
    title(sprintf('Ref Position (Az: %d°, El: %d°)', refAzimuthDeg, refElevationDeg));
    grid on;
    axis equal;

    % Customize tick labels to show degrees instead of radians
    % xticks(deg2rad([-180 -120 -60 0 60 120 180]));
    % xticklabels({'-180°', '-120°', '-60°', '0°', '60°', '120°', '180°'});
    % yticks(deg2rad([-90 -60 -30 0 30 60 90]));
    % yticklabels({'-90°', '-60°', '-30°', '0°', '30°', '60°', '90°'});

    xticks(deg2rad([-180 0 180]));
    xticklabels({' ', ' ', '　'});
    yticks(deg2rad([-90 0 90]));
    yticklabels({' ', ' ', ' '});

    % Plot meridians and parallels
    plotMeridiansAndParallels();

    % Hold off for next plot
    hold off;
end

% Set uniform axis limits
allAxes = findall(gcf, 'Type', 'axes');
for k = 1:numel(allAxes)
    set(allAxes(k), 'XLim', [deg2rad(-30) pi], 'YLim', [-1 1]);
end

% Overall plot title
sgtitle(sprintf('[%s] Mercator Projection Plot of Azimuth and Elevation by Target Position and Stimulus', HRTFtarget), 'Interpreter', 'none');

% Create a common legend in a separate figure
% figure;
% axis off;
legend(plotHandles(1:3), 'Stimulus 0', 'Stimulus 1', 'Reference Point', 'Location', 'best');

    function [x, y] = mercatorProjection(lambda, phi)
        % Mercator projection formula
        % lambda: longitude (radians)
        % phi: latitude (radians)

        x = lambda;
        y = log(tan(pi/4 + phi/2));
    end

    function plotMeridiansAndParallels()
        % Plot meridians
        for lon = -180:30:180
            lonRad = deg2rad(lon);
            latRange = linspace(-85, 85, 100); % Avoiding poles for Mercator
            latRad = deg2rad(latRange);
            [x, y] = mercatorProjection(lonRad * ones(size(latRad)), latRad);
            plot(x, y, 'Color', [0.3 0.3 0.3], 'LineStyle', '--');
            text(x(50), y(50), sprintf('%d°', lon), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', [0.4 0.4 0.4]);
        end

        % Plot parallels
        for lat = -90:30:90
            if lat == 90 || lat == -90
                continue; % Skip poles for Mercator
            end
            latRad = deg2rad(lat);
            lonRange = linspace(-180, 180, 100);
            lonRad = deg2rad(lonRange);
            [x, y] = mercatorProjection(lonRad, latRad * ones(size(lonRad)));
            plot(x, y, 'Color', [0.3 0.3 0.3], 'LineStyle', '--');
            text(x(50), y(50), sprintf('%d°', lat), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', [0.4 0.4 0.4]);
        end
    end
end