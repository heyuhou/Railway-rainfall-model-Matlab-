clear
clc
close all

%% Climate Exposure Mapping

% Read the track data
S = shaperead('Railway_track.shp');

% Extract the coordinates of the midpoint of the track
for i = 1:length(S)
    x = S(i).X(~isnan(S(i).X));
    y = S(i).Y(~isnan(S(i).Y));
    midX(i) = mean(x);
    midY(i) = mean(y);
end

% Read precipitation data
[prGrid, R] = readgeoraster('euro_pr.tif');

% Confirm the spatial range and resolution
disp(R)

% Combine the track data with the climate grid
for i = 1:length(midX)
    [xIntrinsic, yIntrinsic] = geographicToIntrinsic(R, midY(i), midX(i));
    row = round(yIntrinsic);
    col = round(xIntrinsic);

    if row >= 1 && row <= size(prGrid, 1) && col >= 1 && col <= size(prGrid, 2)
        rain_val(i) = prGrid(row, col);
    else
        rain_val(i) = NaN;
    end
end

% Output result
ResultTable = table((1:length(midX))', midX', midY', rain_val', 'VariableNames', {'EdgeID', 'Mid_X', 'Mid_Y', 'Rain'});
writetable(ResultTable, 'rail_with_rain.csv');

% visualization
figure
scatter(midX, midY, 10, rain_val, 'filled')
colormap('hot')
colorbar
cb = colorbar;
cb.Label.String = 'Precipitation (mm/day)'; 
xlabel('Longitude')
ylabel('Latitude')
title('Railway Midpoints with Climate Exposure (maximum one-day precipitation)')
grid on

% Statistics
Max_rain = max(rain_val);
disp(['Maximum: ', num2str(Max_rain), ' mm/day'])
Min_rain = min(rain_val);
disp(['Minimum: ', num2str(Min_rain), ' mm/day'])
Mean_rain = mean(rain_val);
disp(['Mean: ', num2str(Mean_rain), ' mm/day'])
Median_rain = median(rain_val);
disp(['Median: ', num2str(Median_rain), ' mm/day'])
Std_rain = std(rain_val);
disp(['Standard Deviation: ', num2str(Std_rain), ' mm/day'])

% Plot histogram
figure
histogram(rain_val);
xlabel('Precipitation (mm/day)');
ylabel('Number of Railway Segments');
title('Histogram of Precipitation Exposure on Railway Segments');
grid on

%% Fragility Function Development

% Logistic Fragility Function
threshold_log = 40;
k = 0.2197;
failure_prob_log = 1 ./ (1 + exp(-k * (rain_val - threshold_log)));

% Visualization for Logistic
figure
scatter(rain_val, failure_prob_log, 10, 'filled')
xlabel('Precipitation (mm)')
ylabel('Failure Probability')
title('Logistic Fragility Function')
grid on

% % Gaussian Fragility Function
% threshold_gau = 40; % 40mm/day (USGS, 2019; NIC, 2024; UK Government, 2022)
% sigma = 10;
% failure_prob_gau = (1/2) * (1 + (erf((rain_val - threshold_gau)/(sigma * sqrt(2)))));
% 
% % Visualization for Gaussian
% figure
% scatter(rain_val, failure_prob_gau, 10, 'filled')
% xlabel('Precipitation (mm)')
% ylabel('Failure Probability')
% title('Gaussian Fragility Function')
% grid on

%% Failure Scenario Simulation (Monte Carlo) - use logistic

% Parameter Setting
failure_prob_log = failure_prob_log';
N = 1000;
numTracks = length(failure_prob_log);

% Matrix Setting
failure_status = zeros(numTracks, N);

% Monte Carlo Flow
for j = 1:N
    % Generate a random number (0,1)
    randNums = rand(numTracks, 1);

    % Bernoulli criterion, Fail(1), NoFail(0)
    failure_status(:, j) = randNums < failure_prob_log;
end

% Average failure frequency in Monte Carlo
failure_rate_per_track = mean(failure_status, 2);

% Visualization for Monte Carlo
figure
scatter(failure_prob_log, failure_rate_per_track, 'filled')
xlabel('Theoretical failure probability')
ylabel('Simulated average failure frequency (Monte Carlo)')
title('Theoretical failure probability vs Simulated average failure frequency')
grid on
hold on
plot([0,1], [0,1], 'r--') % Ideally
hold off

% Verification of Monte Carlo
residuals = failure_rate_per_track - failure_prob_log;
figure
histogram(residuals, 50, 'Normalization', 'pdf')
xlabel('Residual (Simulated - Theoretical)')
ylabel('Density')
title('Histogram of Simulation Residuals')
grid on

% Heatmap of Simulated Failure Frequencies
figure
scatter(midX, midY, 10, failure_rate_per_track, 'filled')
colormap('hot')
colorbar
cb = colorbar;
cb.Label.String = 'Simulated Failure Frequency'; 
clim([0 1])
xlabel('Longitude')
ylabel('Latitude')
title('Simulated Failure Frequency Heatmap of Railway Segments')
grid on

%% Forecasting Material Needs under Adaptation Scenarios

% Scenario 1: Failure Threshold-Based Allocation
track_length_km = [S.segment_le];
track_length_km = track_length_km(:);
MI = 55;
affected_idx = failure_rate_per_track >= 0.5;
materials_scenario1 = sum(track_length_km(affected_idx) * MI, 'omitnan');
disp(['Total material needed under Scenario 1: ', num2str(materials_scenario1), ' tons']);

% Scenario 2: Probability-Weighted Allocation
materials_scenario2 = sum(track_length_km .* failure_rate_per_track * MI, 'omitnan');
disp(['Total material needed under Scenario 2: ', num2str(materials_scenario2), ' tons']);

% Verification
material_demand = [materials_scenario1, materials_scenario2] / 1e6;
scenarios = {'Scenario 1', 'Scenario 2'};
figure;
bar(material_demand, 'FaceColor', [1 0.6 0]);
xticks(1:2);
xticklabels(scenarios);
ylabel('Additional Stock (Mt)');
title('Additional Stock under Scenario 1 and 2');
grid on;
for i = 1:length(material_demand)
    text(i, material_demand(i) + 0.03, sprintf('%.4f', material_demand(i)), 'HorizontalAlignment', 'center', 'FontSize', 10);
end