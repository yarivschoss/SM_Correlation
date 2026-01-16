function plot_excel_energy_log(xlsxFile)
% plot_excel_energy_log
%   Reads an Excel file with columns:
%       STARTOFINTERVAL_TIMES, EC
%   and plots:
%       X-axis = timestamps
%       Y-axis = energy consumption
%
% Input:
%   xlsxFile - path to Excel file, e.g. "synthetic_customer.xlsx"

    if ~(isstring(xlsxFile) || ischar(xlsxFile))
        error("xlsxFile must be a string or char.");
    end

    % Read Excel file
    T = readtable(xlsxFile);

    % Validate required columns
    if ~all(ismember({'STARTOFINTERVAL_TIMES','EC'}, T.Properties.VariableNames))
        error("Excel file must contain columns: STARTOFINTERVAL_TIMES and EC.");
    end

    % Extract columns
    timeVec = datetime(T.STARTOFINTERVAL_TIMES);
    energy  = T.EC;

    % Plot
    figure('Name','Energy Log from Excel','Color','w');
    plot(timeVec, energy, 'b');
    grid on;
    xlabel("Time (STARTOFINTERVAL_TIMES)");
    ylabel("Energy Consumption (EC)");
    title("Continuous Quarter-Hourly Energy Log");
end
