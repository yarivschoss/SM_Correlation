function plot_excel_energy_log(xlsxFile)
% plot_excel_energy_log
% Reads an Excel file with columns:
% STARTOFINTERVAL_TIMES, EC
% and plots:
% 1) Time series of energy consumption
% 2) Distribution (histogram) and KDE of EC values
%
% Input:
% xlsxFile - path to Excel file, e.g. "synthetic_customer.xlsx"

%% ---------- Input validation ----------
if ~(isstring(xlsxFile) || ischar(xlsxFile))
    error("xlsxFile must be a string or char.");
end

%% ---------- Read Excel ----------
T = readtable(xlsxFile);
if ~all(ismember({'STARTOFINTERVAL_TIMES','EC'}, T.Properties.VariableNames))
    error("Excel file must contain columns: STARTOFINTERVAL_TIMES and EC.");
end

%% ---------- Extract data ----------
timeVec = datetime(T.STARTOFINTERVAL_TIMES);
energy  = T.EC(:);

%% ---------- Statistics ----------
mu  = mean(energy,'omitnan');
med = median(energy,'omitnan');
sd  = std(energy,'omitnan');

%% ---------- Plot ----------
figure('Name','Energy Log Analysis','Color','w');

% ===== Time Series =====
subplot(2,1,1)
plot(timeVec, energy, 'b')
grid on
xlabel("Time")
ylabel("Energy (EC)")
title("Quarter-Hourly Energy Consumption")

% ===== Distribution =====
subplot(2,1,2)

% Changed Normalization to 'pdf' so the histogram and KDE share the same Y-axis scale
histogram(energy, 40, 'Normalization', 'pdf', 'FaceColor',[0.2 0.6 0.8], 'EdgeColor','none', 'DisplayName', 'Histogram');
grid on
xlabel("Energy (EC)")
ylabel("Probability Density") % Updated from Frequency
title("Distribution of Energy Consumption with KDE")
hold on

% ----- Add KDE Curve -----
% Filter out NaNs to prevent ksdensity from failing
valid_energy = energy(~isnan(energy));
[f, xi] = ksdensity(valid_energy);
plot(xi, f, 'Color', [0.85 0.33 0.10], 'LineWidth', 2.5, 'DisplayName', 'KDE Curve');
% -------------------------

xline(mu,  'r','LineWidth',2,'DisplayName','Mean');
xline(med, 'k--','LineWidth',2,'DisplayName','Median');
legend show

% Optional text box with stats
text(0.98,0.95, sprintf('Mean = %.2f\nMedian = %.2f\nStd = %.2f',mu,med,sd), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'BackgroundColor','white','EdgeColor','black');
end