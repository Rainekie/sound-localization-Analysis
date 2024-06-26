function tmpdiffs = diffplot(resArray, referencePoints, HRTFtarget, flag)
if flag(1)
    % Plot the results
    figure;

    % Prepare the tiled layout
    t = tiledlayout(3,2);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    for HRTFs = 1:3
        data = resArray(:,:,HRTFs);
        % Initialize arrays to store differences and standard deviations
        meanDifferences = [];
        stdDifferences = [];
        maeValues = [];
        rmseValues = [];
        targetLabels = cell(size(referencePoints, 1), 1);

        % Loop through each target position
        for i = 1:size(referencePoints, 1)
            currentPosition = referencePoints(i, 1);
            refAzimuth = referencePoints(i, 2);
            refElevation = referencePoints(i, 3);
            refAzimuthRad = deg2rad(refAzimuth);
            refElevationRad = deg2rad(refElevation);

            % Calculate reference position in Cartesian coordinates
            refX = cos(refElevationRad) * cos(refAzimuthRad);
            refY = cos(refElevationRad) * sin(refAzimuthRad);
            refZ = sin(refElevationRad);

            % Filter data for the current target position
            filteredData = data(data(:, 1) == currentPosition, :);
            differences = [];

            % Loop through each stimulus condition
            for j = 0:1
                stimulusData = filteredData(filteredData(:, 2) == j, :);
                azimuthRad = deg2rad(abs(stimulusData(:, 3)));
                elevationRad = deg2rad(stimulusData(:, 4));

                % Calculate user positions in Cartesian coordinates
                x = cos(elevationRad) .* cos(azimuthRad);
                y = cos(elevationRad) .* sin(azimuthRad);
                z = sin(elevationRad);

                % Calculate great-circle distances using Haversine formula
                % Haversine formula parameters
                deltaAzimuth = azimuthRad - refAzimuthRad;
                deltaElevation = elevationRad - refElevationRad;

                % Haversine distance calculation
                a = sin(deltaElevation/2).^2 + cos(refElevationRad) .* cos(elevationRad) .* sin(deltaAzimuth/2).^2;
                c = 2 * atan2(sqrt(a), sqrt(1 - a));
                diff = c; % Since we are on a unit sphere, radius = 1

                differences = [differences; diff];
            end
            % tmpdiffs(:,i, HRTFs) = differences;
            tmpdiffs{i, HRTFs} = differences;

            % Calculate mean and standard deviation of differences
            % コサインとサインを計算
            cosines = cos(differences);
            sines = sin(differences);

            % コサインとサインの平均を計算
            mean_cos = mean(cosines);
            mean_sin = mean(sines);

            % 平均角度を計算
            mean_angle = atan2(mean_sin, mean_cos);
            % 中央値を計算（0から2πの範囲にラッピング）
            wrapped_radians = mod(differences, 2*pi);
            median_angle = median(wrapped_radians);
            
            % ベクトルの長さを計算
            R = sqrt(mean_cos^2 + mean_sin^2);
            
            % 分散と標準偏差を計算
            circular_variance = 1 - R;
            circular_stddev = sqrt(-2 * log(R));

            % Calculate Mean Angular Error (MAE)
            mae = mean(abs(differences));

            % Calculate Root Mean Square Error (RMSE)
            rmse = sqrt(mean(differences.^2));

            meanDifferences = [meanDifferences; mean_angle];
            stdDifferences = [stdDifferences; circular_stddev];
            maeValues = [maeValues; mae];
            rmseValues = [rmseValues; rmse];
            targetLabels{i} = sprintf('Az%d El%d', refAzimuth, refElevation);
        end
        meanDifferences = rad2deg(meanDifferences);
        stdDifferences = rad2deg(stdDifferences);

        % disp('Mean Angular Error (MAE) values:');
        % disp(maeValues);
        % 
        % disp('Root Mean Square Error (RMSE) values:');
        % disp(rmseValues);

        % nexttile;
        % % Bar plot for mean differences
        % b = bar(meanDifferences);
        % hold on;
        % 
        % % Change the color of even-numbered bars
        % barColors = repmat([0 0.4470 0.7410], length(meanDifferences), 1); % Default color
        % barColors(2:2:end, :) = repmat([0.8500 0.3250 0.0980], length(2:2:length(meanDifferences)), 1); % Color for even bars
        % b.FaceColor = 'flat';
        % b.CData = barColors;
        % 
        % % Error bars for standard deviations
        % errorbar(1:length(meanDifferences), meanDifferences, stdDifferences, 'k', 'LineStyle', 'none');
        % 
        % % Set plot labels and title
        % xlabel('Target Position');
        % ylabel('Mean Distance Difference [deg]');
        % ylim([-10, 60])
        % title(sprintf("[%s] each position", HRTFtarget(HRTFs)), 'Interpreter', 'none');
        % set(gca, 'XTick', 1:length(targetLabels), 'XTickLabel', targetLabels);
        % grid on;
        % hold off;
        
        nexttile;
        % Box plot for differences
        data = [];
        group = [];
        for i = 1:length(tmpdiffs)
            data = [data; tmpdiffs{i, HRTFs}];
            group = [group; i * ones(length(tmpdiffs{i, HRTFs}), 1)];
        end
        data = rad2deg(data);

        boxplot(data, group, 'Labels', targetLabels);
        
        % Set plot labels and title
        xlabel('Target Position');
        ylabel('Differences [deg]');
        ylim([-3 75])
        title(sprintf("[%s] each position", HRTFtarget(HRTFs)), 'Interpreter', 'none');
        grid on;

        % mean's mean value
        meanDifferencess = [];
        stdDifferencess = [];
        nums = 1:length(meanDifferences);
        odd_numbers = nums(mod(nums, 2) == 1);
        even_numbers = nums(mod(nums, 2) == 0);
        row = [odd_numbers; even_numbers];
        % コサインとサインを計算
        for j = 1:2
            tmpalldat{j} = cell2mat(tmpdiffs(row(j,:), HRTFs));
            % cosines1 = cos(tmpdiffs(:,row(j,:), HRTFs));
            cosines1 = cos(cell2mat(tmpdiffs(row(j,:), HRTFs)));
            % sines1 = sin(tmpdiffs(:,row(j,:), HRTFs));
            sines1 = sin(cell2mat(tmpdiffs(row(j,:), HRTFs)));

            % コサインとサインの平均を計算
            mean_cos1 = mean(cosines1, 'all');
            mean_sin1 = mean(sines1, 'all');

            % 平均角度を計算
            mean_angle1 = atan2(mean_sin1, mean_cos1);

            % 中央値を計算（0から2πの範囲にラッピング）
            % wrapped_radians1 = mod(tmpdiffs(:,row(j,:), HRTFs), 2*pi);
            wrapped_radians1 = mod(cell2mat(tmpdiffs(row(j,:), HRTFs)), 2*pi);
            median_angle1 = median(wrapped_radians1);

            % ベクトルの長さを計算
            R1 = sqrt(mean_cos1^2 + mean_sin1^2);

            % 分散と標準偏差を計算
            circular_variance1 = 1 - R1;
            circular_stddev1 = sqrt(-2 * log(R1));

            meanDifferencess = [meanDifferencess; mean_angle1];
            stdDifferencess = [stdDifferencess; circular_stddev1];
        end
        meanDifferencess = rad2deg(meanDifferencess);
        stdDifferencess = rad2deg(stdDifferencess);

        % nexttile;
        % % Bar plot for mean differences
        % b = bar(meanDifferencess);
        % hold on;
        % 
        % % Change the color of even-numbered bars
        % barColors = repmat([0 0.4470 0.7410], length(meanDifferencess), 1); % Default color
        % barColors(2:2:end, :) = repmat([0.8500 0.3250 0.0980], length(2:2:length(meanDifferencess)), 1); % Color for even bars
        % b.FaceColor = 'flat';
        % b.CData = barColors;
        % 
        % % Error bars for standard deviations
        % errorbar(1:length(meanDifferencess), meanDifferencess, stdDifferencess, 'k', 'LineStyle', 'none');
        % 
        % % Set plot labels and title
        % xlabel('Target Position');
        % ylabel('Mean Distance Difference [deg]');
        % ylim([-10, 40])
        % title(sprintf("[%s] vertical comparison", HRTFtarget(HRTFs)), 'Interpreter', 'none');
        % set(gca, 'XTick', 1:2, 'XTickLabel', {'Horizontal', 'Upper'});
        % grid on;
        % hold off;

        % ボックスプロット表示用のデータとグループインデックスを準備
        data = [];
        group = [];
        for i = 1:length(tmpalldat)
            data = [data; tmpalldat{i}];
            group = [group; i * ones(length(tmpalldat{i}), 1)];
        end
        data = rad2deg(data);

        % ボックスプロットを作成
        nexttile;
        boxplot(data, group);
        xlabel('Target Position');
        ylabel('Mean Distance Difference [deg]');
        ylim([-3 75])
        title(sprintf("[%s] vertical comparison", HRTFtarget(HRTFs)), 'Interpreter', 'none');
        set(gca, 'XTick', 1:2, 'XTickLabel', {'Horizontal', 'Upper'});
        grid on;
    end

    title(t, 'Mean Great-Circular Distance Difference from Reference Positions with Variability');
end

if flag(2)
    % Plot the results
    figure;

    % Prepare the tiled layout
    t = tiledlayout(3,2);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    for HRTFs = 1:3
        data = resArray(:,:,HRTFs);
        % Initialize arrays to store differences and standard deviations
        meanAzimuthDifferences = [];
        stdAzimuthDifferences = [];
        meanElevationDifferences = [];
        stdElevationDifferences = [];
        targetLabels = cell(size(referencePoints, 1), 1);

        % Loop through each target position
        for i = 1:size(referencePoints, 1)
            currentPosition = referencePoints(i, 1);
            refAzimuth = referencePoints(i, 2);
            refElevation = referencePoints(i, 3);
            refAzimuthRad = deg2rad(refAzimuth);
            refElevationRad = deg2rad(refElevation);

            % Calculate reference position in Cartesian coordinates
            refX = cos(refElevationRad) * cos(refAzimuthRad);
            refY = cos(refElevationRad) * sin(refAzimuthRad);
            refZ = sin(refElevationRad);

            % Filter data for the current target position
            filteredData = data(data(:, 1) == currentPosition, :);
            azimuthDifferences = [];
            elevationDifferences = [];

            % Loop through each stimulus condition
            for j = 0:1
                stimulusData = filteredData(filteredData(:, 2) == j, :);
                azimuth = abs(stimulusData(:, 3)); % Azimuth in degrees
                elevation = stimulusData(:, 4); % Elevation in degrees

                % Calculate differences
                azimuthDiff = abs(azimuth - refAzimuth);
                elevationDiff = abs(elevation - refElevation);
                azimuthDifferences = [azimuthDifferences; azimuthDiff];
                elevationDifferences = [elevationDifferences; elevationDiff];
            end

            % Calculate mean and standard deviation of differences
            meanAzimuthDiff = mean(azimuthDifferences);
            stdAzimuthDiff = std(azimuthDifferences);
            meanElevationDiff = mean(elevationDifferences);
            stdElevationDiff = std(elevationDifferences);
            meanAzimuthDifferences = [meanAzimuthDifferences; meanAzimuthDiff];
            stdAzimuthDifferences = [stdAzimuthDifferences; stdAzimuthDiff];
            meanElevationDifferences = [meanElevationDifferences; meanElevationDiff];
            stdElevationDifferences = [stdElevationDifferences; stdElevationDiff];
            targetLabels{i} = sprintf('Az%d El%d', refAzimuth, refElevation);
        end

        % Plot for Azimuth Differences
        nexttile;
        bAz = bar(meanAzimuthDifferences);
        hold on;

        % Change the color of even-numbered bars
        barColorsAz = repmat([0 0.4470 0.7410], length(meanAzimuthDifferences), 1); % Default color
        barColorsAz(2:2:end, :) = repmat([0.8500 0.3250 0.0980], length(2:2:length(meanAzimuthDifferences)), 1); % Color for even bars
        bAz.FaceColor = 'flat';
        bAz.CData = barColorsAz;

        % Error bars for standard deviations
        errorbar(1:length(meanAzimuthDifferences), meanAzimuthDifferences, stdAzimuthDifferences, 'k', 'LineStyle', 'none');

        % Set plot labels and title
        xlabel('Target Position');
        ylabel('Mean Azimuth Difference (degrees)');
        ylim([-5, 30])
        title(sprintf("[%s] Azimuth Differences", HRTFtarget(HRTFs)), 'Interpreter', 'none');
        set(gca, 'XTick', 1:length(targetLabels), 'XTickLabel', targetLabels);
        grid on;
        hold off;

        % Plot for Elevation Differences
        nexttile;
        bEl = bar(meanElevationDifferences);
        hold on;

        % Change the color of even-numbered bars
        barColorsEl = repmat([0 0.4470 0.7410], length(meanElevationDifferences), 1); % Default color
        barColorsEl(2:2:end, :) = repmat([0.8500 0.3250 0.0980], length(2:2:length(meanElevationDifferences)), 1); % Color for even bars
        bEl.FaceColor = 'flat';
        bEl.CData = barColorsEl;

        % Error bars for standard deviations
        errorbar(1:length(meanElevationDifferences), meanElevationDifferences, stdElevationDifferences, 'k', 'LineStyle', 'none');

        % Set plot labels and title
        xlabel('Target Position');
        ylabel('Mean Elevation Difference (degrees)');
        ylim([-5, 30])
        title(sprintf("[%s] Elevation Differences", HRTFtarget(HRTFs)), 'Interpreter', 'none');
        set(gca, 'XTick', 1:length(targetLabels), 'XTickLabel', targetLabels);
        grid on;
        hold off;
    end

    title(t, 'Mean Azimuth and Elevation Differences from Reference Positions with Variability');
end

if flag(3)
    % Plot the results
    figure;

    % Prepare the tiled layout
    t = tiledlayout(3,1);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    for HRTFs = 1:3
        data = resArray(:,:,HRTFs);
        % Initialize arrays to store mean and standard deviations of reaction times
        meanReactionTimes = [];
        stdReactionTimes = [];
        targetLabels = cell(size(referencePoints, 1), 1);

        % Loop through each target position
        for i = 1:size(referencePoints, 1)
            currentPosition = referencePoints(i, 1);

            % Filter data for the current target position
            filteredData = data(data(:, 1) == currentPosition, :);
            reactionTimes = [];

            % Loop through each stimulus condition
            for j = 0:1
                stimulusData = filteredData(filteredData(:, 2) == j, :);
                reactionTimes = [reactionTimes; stimulusData(:, end)]; % Reaction time is the last column
            end

            % Calculate mean and standard deviation of reaction times
            meanRT = mean(reactionTimes);
            stdRT = std(reactionTimes);
            meanReactionTimes = [meanReactionTimes; meanRT];
            stdReactionTimes = [stdReactionTimes; stdRT];
            targetLabels{i} = sprintf('Az%d El%d', referencePoints(i, 2), referencePoints(i, 3));
        end
        clear data

        % Bar plot for mean reaction times
        nexttile;
        bRT = bar(meanReactionTimes);
        hold on;

        % Change the color of even-numbered bars
        barColorsRT = repmat([0 0.4470 0.7410], length(meanReactionTimes), 1); % Default color
        barColorsRT(2:2:end, :) = repmat([0.8500 0.3250 0.0980], length(2:2:length(meanReactionTimes)), 1); % Color for even bars
        bRT.FaceColor = 'flat';
        bRT.CData = barColorsRT;

        % Error bars for standard deviations
        errorbar(1:length(meanReactionTimes), meanReactionTimes, stdReactionTimes, 'k', 'LineStyle', 'none');

        % Set plot labels and title
        xlabel('Target Position');
        ylabel('Mean Reaction Time (s)');
        ylim([-5 20])
        title(sprintf("[%s]", HRTFtarget(HRTFs)), 'Interpreter', 'none');
        set(gca, 'XTick', 1:length(targetLabels), 'XTickLabel', targetLabels);
        grid on;
        hold off;
    end

    title(t, 'Mean Reaction Time for Different Conditions with Variability');
end
end