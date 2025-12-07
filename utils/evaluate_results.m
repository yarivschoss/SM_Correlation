function evaluate_results(results, cfg)
    % evaluate_results – print summary, compute metrics and generate plots
    %
    % Examples:
    %   - Size and sparsity of assignmentMatrix
    %   - Accuracy metrics compared to ground truth (if available)
    %   - Plots of costHistory / residualHistory

    disp("=== Energy Balance Result ===");
    disp(size(results.energyBalance.assignmentMatrix));

    % TODO: add accuracy calculation, plotting of costHistory, etc.
end
