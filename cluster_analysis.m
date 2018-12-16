function cluster_analysis(write_pdf)
if nargin == 0
    write_pdf = 0;
end
data = read_all();

Y = 2007: 2012;
K = 2: 1: 10; % # of clusters
% f_cluster = 'kmedoids';
f_cluster = 'kmeans';
s_cluster   = nan(length(K), length(Y)); % Silhouette metric
d_cluster   = nan(length(K), length(Y)); % Dunn metric
idxC = nan(data.nI, length(K), size(Y, 1));
for j = 1: length(Y)
    y = Y(j);
    xa_t = data.va(data.yyyy==y, :)'; % Let's try wind speed
    for i = 1: length(K)
        k = K(i);
        [idxC(:, i, j), ~] = feval(...
            f_cluster, xa_t, k,...
            'Replicates',10,...
            'Distance', 'correlation');
        [s, ~] = silhouette(xa_t, idxC(:, i, j), 'correlation');
        s_cluster(i, j) = mean(s);
        distM=squareform(pdist(xa_t));
        d_cluster(i, j) = dunns(k, distM, idxC(:, i, j));
    end
    fig_map(idxC(:, K==2));
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
end

figure();
h = plot(K, s_cluster, '-s');
set(h(1), 'LineStyle','-',  'LineWidth', 2, 'Color', 'k');
set(h(2), 'LineStyle','--', 'LineWidth', 2, 'Color', 'k');
set(h(3), 'LineStyle','-.', 'LineWidth', 2, 'Color', 'k');
set(h(4), 'LineStyle','-',  'LineWidth', 2, 'Color', 'b');
set(h(5), 'LineStyle','--', 'LineWidth', 2, 'Color', 'b');
set(h(6), 'LineStyle','-.', 'LineWidth', 2, 'Color', 'b');
Ystr{j} = int2str(y);
for j = 1: length(Y)
    y = Y(j);
    Ystr{j} = int2str(y);
end
legend(h, Ystr);
xlabel('# of clusters');
set(findall(gcf,'-property','FontSize'),'FontSize',24);
if write_pdf == 1
    h = gcf();
    set(h, 'Units', 'Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,'Silhouette','-dpdf','-r0');
end
title('Silhouette');

figure();
plot(K, d_cluster, '-s');
xlabel('# of clusters');
title('Dunn');

figure();
set(gcf, 'Position', [0, 0, 500, 600]);
fig.PaperPositionMode = 'auto';
for j = 1:length(Y)
    y = Y(j);
    r_va = corrcoef(data.va(data.yyyy==y, :));
    ax = subplot(3, 2, j);
    pos0 = ax.Position;
    if mod(j, 2) == 0
        ax.Position = [pos0(1)-0.05, pos0(2), pos0(3)*0.95, pos0(4)];
    else
        ax.Position = [pos0(1)-0.03, pos0(2), pos0(3)*0.95, pos0(4)];
    end
    imagesc(r_va);
    colormap jet;
    title(y);
end
pos0 = ax.Position;
colorbar('Position', [pos0(1)+pos0(3)+0.05  pos0(2)  0.03  pos0(4)*3.8]);

figure();
r_va = corrcoef(data.va);
imagesc(r_va);
colormap jet;
colorbar;
set(findall(gcf,'-property','FontSize'),'FontSize',24);
if write_pdf == 1
    h = gcf();
    set(h, 'Units', 'Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,'corr','-dpdf','-r0');
end

figure();
idx2 = squeeze(idxC(:, 1, :));
% idx2 = idx2 - 1;
% idx2(:, 1) = -1.*(idx2(:, 1) - 1);
% idx2(:, 5) = -1.*(idx2(:, 5) - 1);
[r, c] = size(idx2');
imagesc((1:c)+0.5, (1:r)+0.5, idx2');
colormatrix = [255, 137, 137; 138, 199, 252]./255;
colormap(colormatrix);
ax = gca;
set(ax, 'XTick', 1:(c+1), 'YTick', 1:(r+1),...
    'XLim', [1 c+1], 'YLim', [1 r+1],...
    'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
set(ax, 'YDir', 'normal');
set(findall(gcf,'-property','FontSize'),'FontSize',14);
for i = 1: length(ax.XTickLabel)
    if mod(i, 5) == 0
        ax.XTickLabel{i} = int2str(i);
    else
        ax.XTickLabel{i} = '';
    end
end
for i = 1: length(Y)
    ax.YTickLabel{i} = int2str(Y(i));
end
ax.YTickLabel{end} = '';
if write_pdf == 1
    h = gcf();
    set(h, 'Units', 'Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,'2cluster','-dpdf','-r0');
end
end