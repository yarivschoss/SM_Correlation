function fig = plot_topology_gt_vs_method(results, varargin)
% plot_topology_gt_vs_method
% Visual comparison: GT vs predicted topology with error metrics and mismatch highlighting.

% ---------- Parse Name-Value optional inputs ----------
p = inputParser;
p.addParameter('Title', "", @(x)isstring(x)||ischar(x));
p.addParameter('GridSpacing', 2.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('SquareSize', 1.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('FontSize', 9, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('DrawLinks', false, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('ShowLegend', true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

% ---------- Extract data from results struct ----------
if iscell(results.algorithm)
    methodName = string(results.algorithm{1});
else
    methodName = string(results.algorithm);
end

customersDir = string(results.customersDir);
isSon_pred   = logical(results.isSon);

% ---------- Get Ground Truth (GT) ----------
if isfield(results.classification, 'gt_isSon') && ~isempty(results.classification.gt_isSon)
    isSon_gt = logical(results.classification.gt_isSon);
else
    metaPath = fullfile(customersDir, "customers_metadata.xlsx");
    if ~isfile(metaPath), error("Ground Truth not found in results or Excel."); end
    M = readtable(metaPath);
    vars = lower(string(M.Properties.VariableNames));
    idx_is_son = find(vars=="is_son" | vars=="isson", 1);
    isSon_gt = logical(M{:, idx_is_son});
end

isSon_gt = isSon_gt(:);
isSon_pred = isSon_pred(:);
N = numel(isSon_gt);

% ---------- Calculate Error Metrics ----------
mismatches = (isSon_gt ~= isSon_pred); % Logical vector of errors
numErrors = sum(mismatches);
errorRate = (numErrors / N) * 100;

% ---------- Figure Setup ----------
fig = figure('Color','w','Name',"Topology Comparison - " + methodName, 'Position', [100 100 1100 550]);
xy = buildGridWithHole(N, opt.GridSpacing);

% Colors
cSon    = [1.00 0.55 0.10];  % orange
cOrphan = [0.20 0.50 0.95];  % blue
cTx     = [0.20 0.20 0.20];  % transformer

% Title with stats
mainTitle = opt.Title; if strlength(mainTitle)==0, mainTitle = "Topology: GT vs " + methodName; end
fullTitle = sprintf('%s\nMismatches: %d / %d (%.1f%% Error Rate)', mainTitle, numErrors, N, errorRate);
sgtitle(fullTitle, 'FontWeight','bold', 'FontSize', 13);

% Panel 1: Ground Truth
ax1 = subplot(1,2,1); hold(ax1,'on'); axis(ax1,'equal'); axis(ax1,'off');
title(ax1, "Ground Truth (Actual System)");
draw_panel_grid(ax1, xy, isSon_gt, false(N,1), opt.SquareSize, opt.FontSize, cSon, cOrphan, cTx, opt.DrawLinks);

% Panel 2: Prediction with Error Highlighting
ax2 = subplot(1,2,2); hold(ax2,'on'); axis(ax2,'equal'); axis(ax2,'off');
title(ax2, methodName + " (Algorithm Prediction)");
draw_panel_grid(ax2, xy, isSon_pred, mismatches, opt.SquareSize, opt.FontSize, cSon, cOrphan, cTx, opt.DrawLinks);

% Legend
if opt.ShowLegend
    axes(ax2); %#ok<LAXES>
    h1 = plot(nan,nan,'s','MarkerFaceColor',cSon,'MarkerEdgeColor','k','MarkerSize',10);
    h2 = plot(nan,nan,'s','MarkerFaceColor',cOrphan,'MarkerEdgeColor','k','MarkerSize',10);
    h_err = plot(nan,nan,'s','MarkerFaceColor','none','MarkerEdgeColor','r','LineWidth',2,'MarkerSize',12);
    legend([h1 h2 h_err], {'Son','Orphan','Mismatch (Error)'}, 'Location','southoutside','Orientation','horizontal');
end
end

% ==================== INTERNAL HELPER FUNCTIONS ====================


function xy = buildGridWithHole(N, spacing)
    nCols = ceil(sqrt(N+1));
    nRows = ceil((N+1)/nCols);
    [xg, yg] = meshgrid(1:nCols, 1:nRows);
    xy_all = [xg(:), yg(:)];
    xy_all(:,1) = xy_all(:,1) - mean(xy_all(:,1));
    xy_all(:,2) = xy_all(:,2) - mean(xy_all(:,2));
    [~, idxCenter] = min(sum(xy_all.^2,2));
    xy_all(idxCenter,:) = [];
    xy = xy_all(1:N,:) * spacing;
end

function draw_panel_grid(ax, xy, isSon, mismatches, sq, fs, cSon, cOrphan, cTx, drawLinks)
    N = numel(isSon);
    if drawLinks
        for i = 1:N
            plot(ax, [0 xy(i,1)], [0 xy(i,2)], '-', 'Color',[0.9 0.9 0.9], 'LineWidth',0.5);
        end
    end
    % Transformer
    plot(ax, 0,0,'o','MarkerSize',22,'MarkerFaceColor',cTx,'MarkerEdgeColor','k');
    text(ax, 0, -1, "TX", 'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
    
    half = sq/2;
    for i = 1:N
        x = xy(i,1); y = xy(i,2);
        fc = cOrphan; if isSon(i), fc = cSon; end
        
        % Highlight errors with a thick red border
        edgeCol = 'k'; lw = 0.5;
        if mismatches(i), edgeCol = 'r'; lw = 2.5; end
        
        patch(ax, [x-half x+half x+half x-half], [y-half y-half y+half y+half], ...
              fc, 'EdgeColor', edgeCol, 'LineWidth', lw);
          
        tc = [1 1 1]; if isSon(i), tc = [0 0 0]; end
        text(ax, x, y, num2str(i), 'HorizontalAlignment','center', 'Color', tc, 'FontWeight','bold', 'FontSize', fs);
    end
    pad = 2;
    xlim(ax, [min(xy(:,1))-pad, max(xy(:,1))+pad]);
    ylim(ax, [min(xy(:,2))-pad, max(xy(:,2))+pad]);
end