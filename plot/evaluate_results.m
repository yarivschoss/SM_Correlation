function evaluate_results(results)
    % 1. Extract and sanitize the algorithm name
    if iscell(results.algorithm)
        algoName = char(results.algorithm{1});
    else
        algoName = char(results.algorithm);
    end
    
    cm = results.classification.cm;
    classif = results.classification;
    tau_chosen = classif.tau;
    
    % Prepare Confusion Matrix Data
    cm_matrix = [cm.TN, cm.FP; cm.FN, cm.TP];
    
    % Create Figure
    figure('Name', ['Results - ' algoName], 'Color', 'w', 'Position', [100, 100, 1100, 500]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact');
    title(t, ['Performance Analysis: ' algoName], 'FontSize', 16, 'FontWeight', 'bold');

    %% --- Tile 1: Confusion Matrix ---
    nexttile;
    h = confusionchart(cm_matrix, {'Negative', 'Positive'}, 'Title', 'Confusion Matrix');
    
    % Display metrics below the chart
    metricText = sprintf('Precision: %.2f | Recall: %.2f | F1: %.2f', ...
        classif.Precision, classif.Recall, classif.F1);
    xlabel(h, metricText); 

    %% --- Tile 2: ROC Curve ---
    nexttile;
    if isfield(classif, 'ROC_TPR') && isfield(classif, 'ROC_FPR')
        TPR = classif.ROC_TPR;
        FPR = classif.ROC_FPR;
        
        % Sort vectors to ensure correct AUC calculation
        pairs = unique([FPR(:) TPR(:)], 'rows');
        pairs = sortrows(pairs,1);

        FPR_sorted = pairs(:,1);
        TPR_sorted = pairs(:,2);
        
        % 1. Calculate AUC using trapezoidal numerical integration
        auc_val = trapz(FPR_sorted, TPR_sorted);
    else
        % Fallback for single point
        auc_val = NaN;
        TPR = classif.Recall; 
        FPR = cm.FP / max(1, cm.FP + cm.TN);
        FPR_sorted = [0; FPR; 1];
        TPR_sorted = [0; TPR; 1];
    end
    
    % Plot the ROC line
    plot(FPR_sorted, TPR_sorted, 'LineWidth', 2.5, 'Color', [0 0.45 0.74]);
    hold on;
    plot([0 1], [0 1], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2); % Diagonal
    
    % 2. Mark the CHOSEN tau on the curve
    % The point corresponding to our selected tau from the confusion matrix
    pt_FPR = cm.FP / max(1, cm.FP + cm.TN);
    pt_TPR = classif.Recall; % This is the TPR at the chosen tau
    
    plot(pt_FPR, pt_TPR, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
    
    % Add annotation for the chosen tau
    text(pt_FPR + 0.03, pt_TPR - 0.05, sprintf('\\tau_{chosen} = %.2f', tau_chosen), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', 'r');
    
    % Formatting
    grid on; axis square;
    xlim([-0.01 1.01]); ylim([-0.01 1.01]);
    xlabel('False Positive Rate (FPR)');
    ylabel('True Positive Rate (TPR)');
    
    % Display AUC in the title or legend
    title(sprintf('ROC Curve (AUC = %.3f)', auc_val));
    legend('Algorithm ROC', 'Random Guess', 'Chosen Threshold', 'Location', 'southeast');
    
    hold off;
end