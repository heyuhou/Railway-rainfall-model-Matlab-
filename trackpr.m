clear
clc
close all

%% Climate Exposure Mapping

% Read the track data
S = shaperead('track.shp');

% Extract the coordinates of the midpoint of the track
for i = 1:length(S)
    x = S(i).X(~isnan(S(i).X));
    y = S(i).Y(~isnan(S(i).Y));
    midX(i) = mean(x);
    midY(i) = mean(y);
end

% Read precipitation data
[prGrid, R] = readgeoraster('pr.tif');

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
colormap('jet')
colorbar
xlabel('Longitude')
ylabel('Latitude')
title('Railway Midpoints with Climate Exposure')
grid on

%% Fragility Function Development

% Logistic Fragility Function
threshold_log = 40; % 40mm/day (USGS, 2019; NIC, 2024; UK Government, 2022)
k = 0.5;
failure_prob_log = 1 ./ (1 + exp(-k * (rain_val - threshold_log)));

% Visualization for Logistic
figure
scatter(rain_val, failure_prob_log, 10, 'filled')
xlabel('Precipitation (mm)')
ylabel('Failure Probability')
title('Logistic Fragility Function')
grid on

% Gaussian Fragility Function
threshold_gau = 40; % 40mm/day (USGS, 2019; NIC, 2024; UK Government, 2022)
sigma = 5;
failure_prob_gau = (1/2) * (1 + (erf((rain_val - threshold_gau)/(sigma * sqrt(2)))));

% Visualization for Gaussian
figure
scatter(rain_val, failure_prob_gau, 10, 'filled')
xlabel('Precipitation (mm)')
ylabel('Failure Probability')
title('Gaussian Fragility Function')
grid on

%% Failure Scenario Simulation (Monte Carlo) - use logistic

% Parameter Setting
failure_prob_log = failure_prob_log';
N = 500;
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