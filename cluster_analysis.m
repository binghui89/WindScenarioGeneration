function cluster_analysis()
dirhome = pwd;
dirwork = 'texas2000_wind_plus';
%% Read metafile
metafile = 'C:\Users\bxl180002\OneDrive\Tmp_RampWind\code_before_github\WIND\wtk_site_metadata.csv';
fid = fopen(metafile);
str_format = '%s %s %s %s %s %s %s %s %s %s %s %s';
textscan(fid, str_format, 1, 'Delimiter', ',');
str_format = '%d %f %f %s %s %f %s %f %f %f %f %s';
C = textscan(fid, str_format, 'Delimiter', ',', 'EmptyValue',nan);
fclose(fid);

%% Read site ID and find out site capacity
filenames = dir(dirwork);
nI = 0;
for i = 1: length(filenames)
    if filenames(i).isdir
        continue;
    else
        nI = nI + 1;
    end
end
sid = nan(length(filenames), 1);
for i = 1: length(filenames)
    if filenames(i).isdir
        sid(i) = nan;
    else
        [~, name, ~] = fileparts(filenames(i).name);
        sid(i) = str2num(name);
    end
end
sid(isnan(sid)) = [];
sid = sort(sid); % sid in ascending order
cap = nan(size(sid));
for i = 1: length(sid)
    cap(i) = C{8}(C{1}==sid(i));
end

%% Read wind power
cd(dirwork);
fname = strcat(int2str(sid(i)), '.csv');
M = csvread(fname, 1, 0);
yyyy = M(:, 1);
mm   = M(:, 2);
dd   = M(:, 3);
hh   = M(:, 4);
mi   = M(:, 5); 

nT = size(M, 1);
nI = length(sid);

xa_MW = nan(nT, nI); % Actual wind power in MW
xf_MW = nan(nT, nI); % Forecasted wind power in MW
xa    = nan(nT, nI); % Actual wind power in p.u.
xf    = nan(nT, nI); % Forecasted wind power in p.u.
va    = nan(nT, nI); % Actual wind speed in m/s
for i = 1: nI
    fname = strcat(int2str(sid(i)), '.csv');
    M = csvread(fname, 1, 0);
    xa_MW(:, i) = M(1:nT, 8);
    xf_MW(:, i) = M(1:nT, 11);
    xa(:, i)    = M(1:nT, 8)./cap(i);
    xf(:, i)    = M(1:nT, 11)./cap(i);
    va(:, i)    = M(1:nT, 10);
end
cd(dirhome);

%% Clustering analysis
% Y = 2007: 2012;
Y = 2007;
K = 2: 2: 40; % # of clusters
% f_cluster = 'kmedoids';
f_cluster = 'kmeans';
s_cluster   = nan(length(K), length(Y)); % Silhouette metric
d_cluster   = nan(length(K), length(Y)); % Dunn metric
i_selected = (yyyy==2012)&(mm==5);
for j = 1: length(Y)
    y = Y(j);
%     xa_t = xa(yyyy==y, :)';
%     xa_t = va(yyyy==y, :)'; % Let's try wind speed
    xa_t = va(:, :)';
    idx_kmeans = nan(size(xa_t, 1), length(K));
    for i = 1: length(K)
        k = K(i);
        [idx_kmeans(:, i), ~] = feval(...
            f_cluster, xa_t, k,...
            'Replicates',10,...
            'Distance', 'correlation');
        [s, ~] = silhouette(xa_t, idx_kmeans(:, i), 'correlation');
        s_cluster(i, j) = mean(s);
        distM=squareform(pdist(xa_t));
        d_cluster(i, j) = dunns(k, distM, idx_kmeans(:, i));
    end
    fig_map(idx_kmeans(:, K==2));
end
figure();
plot(K, s_cluster, '-s');
title('Silhouette');
xlabel('# of clusters');
figure();
plot(K, d_cluster, '-s');
xlabel('# of clusters');
title('Dunn');
end