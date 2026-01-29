function fig = plot_topology_gt_vs_method(isSon_pred, customersDir, methodName, varargin)
% plot_topology_gt_vs_method
% Visual comparison: GT vs predicted topology (sons/orphans) around transformer (TX).
% Layout: Square grid with a center hole for TX. Each customer is a colored square + index label.
%
% Inputs:
%   isSon_pred   : logical/0-1 vector (Nx1) predicted sons (true=son, false=orphan)
%   customersDir : folder containing customers_metadata.xlsx (with is_son or isSon column)
%   methodName   : string like "EnergyBalance", "Optimization", "ECPC"
%
% Name-Value (optional):
%   'Title'        : custom figure title
%   'GridSpacing'  : distance between meters in grid units (default 2.0)
%   'SquareSize'   : square edge length in axis units (default 1.5)
%   'FontSize'     : label font size (default 9)
%   'DrawLinks'    : draw faint lines from TX to meters (default false)
%   'ShowLegend'   : show legend (default true)

p = inputParser;
p.addParameter('Title', "", @(x)isstring(x)||ischar(x));
p.addParameter('GridSpacing', 2.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('SquareSize', 1.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('FontSize', 9, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('DrawLinks', false, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('ShowLegend', true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

customersDir = string(customersDir);
methodName   = string(methodName);

% ---------- Load GT is_son from metadata ----------
metaPath = fullfile(customersDir, "customers_metadata.xlsx");
assert(isfile(metaPath), "Missing customers_metadata.xlsx in: %s", customersDir);

M = readtable(metaPath);
vars = lower(string(M.Properties.VariableNames));
idx_is_son = find(vars=="is_son" | vars=="isson", 1);
assert(~isempty(idx_is_son), "customers_metadata.xlsx must contain column is_son (or isSon).");

isSon_gt = logical(M{:, idx_is_son});
isSon_gt = isSon_gt(:);

% ---------- Validate sizes ----------
isSon_pred = logical(isSon_pred(:));
N = numel(isSon_gt);
assert(numel(isSon_pred)==N, "Length mismatch: isSon_pred=%d, GT=%d", numel(isSon_pred), N);

% ---------- Build grid positions with center hole for TX ----------
xy = buildGridWithHole(N, opt.GridSpacing);

% ---------- Colors ----------
cSon    = [1.00 0.55 0.10];  % orange
cOrphan = [0.20 0.50 0.95];  % blue
cTx     = [0.20 0.20 0.20];  % transformer (dark)

% ---------- Figure ----------
fig = figure('Color','w','Name',"Topology: GT vs "+methodName);
ttl = string(opt.Title);
if strlength(ttl)==0
    ttl = "GT vs " + methodName;
end
sgtitle(ttl, 'FontWeight','bold');

ax1 = subplot(1,2,1); hold(ax1,'on'); axis(ax1,'equal'); axis(ax1,'off');
title(ax1, "Ground Truth (GT)");

ax2 = subplot(1,2,2); hold(ax2,'on'); axis(ax2,'equal'); axis(ax2,'off');
title(ax2, methodName + " (Pred)");

% ---------- Draw both panels ----------
draw_panel_grid(ax1, xy, isSon_gt,  opt.SquareSize, opt.FontSize, cSon, cOrphan, cTx, opt.DrawLinks);
draw_panel_grid(ax2, xy, isSon_pred,opt.SquareSize, opt.FontSize, cSon, cOrphan, cTx, opt.DrawLinks);

% ---------- Legend ----------
if opt.ShowLegend
    axes(ax2); %#ok<LAXES>
    h1 = plot(nan,nan,'s','MarkerFaceColor',cSon,'MarkerEdgeColor','k','MarkerSize',10);
    h2 = plot(nan,nan,'s','MarkerFaceColor',cOrphan,'MarkerEdgeColor','k','MarkerSize',10);
    h3 = plot(nan,nan,'o','MarkerFaceColor',cTx,'MarkerEdgeColor','k','MarkerSize',12);
    legend([h1 h2 h3], {'son','orphan','transformer'}, ...
        'Location','southoutside','Orientation','horizontal');
end

end

% ==================== helpers ====================

function xy = buildGridWithHole(N, spacing)
% Build a near-square grid with one removed point closest to (0,0) for TX.
nCols = ceil(sqrt(N+1));     % +1 because we remove one (center) later
nRows = ceil((N+1)/nCols);

[xg, yg] = meshgrid(1:nCols, 1:nRows);
xy_all = [xg(:), yg(:)];

% center around (0,0)
xy_all(:,1) = xy_all(:,1) - mean(xy_all(:,1));
xy_all(:,2) = xy_all(:,2) - mean(xy_all(:,2));

% remove the point closest to (0,0) => hole for TX
[~, idxCenter] = min(sum(xy_all.^2,2));
xy_all(idxCenter,:) = [];

% take first N
xy = xy_all(1:N,:) * spacing;

% (Optional) sort by angle or keep natural order.
% Keeping natural order means indices will fill row-by-row; it's usually clearer.
end

function draw_panel_grid(ax, xy, isSon, sq, fs, cSon, cOrphan, cTx, drawLinks)
N = numel(isSon);

% optional faint links
if drawLinks
    for i = 1:N
        plot(ax, [0 xy(i,1)], [0 xy(i,2)], '-', 'Color',[0.88 0.88 0.88], 'LineWidth',0.5);
    end
end

% draw TX at center
plot(ax, 0,0,'o','MarkerSize',20,'MarkerFaceColor',cTx,'MarkerEdgeColor','k','LineWidth',1.0);
text(ax, 0, -0.85, "TX", 'HorizontalAlignment','center','FontWeight','bold','Color',cTx);

% draw squares + indices
half = sq/2;
for i = 1:N
    x = xy(i,1); y = xy(i,2);
    fc = cOrphan;
    if isSon(i), fc = cSon; end

    Xp = [x-half x+half x+half x-half];
    Yp = [y-half y-half y+half y+half];
    patch(ax, Xp, Yp, fc, 'EdgeColor','k', 'LineWidth',0.9);

    % choose text color (white on blue, black on orange)
    if isSon(i)
        tc = [0 0 0];
    else
        tc = [1 1 1];
    end

    text(ax, x, y, num2str(i), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize', fs, ...
        'Color', tc, ...
        'FontWeight','bold');
end

% limits
pad = 1.5;
xlim(ax, [min(xy(:,1))-pad, max(xy(:,1))+pad]);
ylim(ax, [min(xy(:,2))-pad, max(xy(:,2))+pad]);
end
