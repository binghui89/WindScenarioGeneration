% This script reads wind power data of all 45 sites, note the time interval
% includes everything from 2007 to 2012.
% nI: # of sites, nT: # of hours
% sid: site ID in ascending order, cap: capacity in the order of sid
% yyyy: year, mm: month, dd: day, hh: hour, mi: minute
% va: actual wind speed in m/s
% xa_MW: actual site wind power in MW, xa: actual site wind power in p.u.
% xf_MW: forecasted site wind power in MW
% xf: forecasted site wind power in p.u.
function data = read_all()
dirhome = pwd;
dirwork = 'texas2000_wind_plus';
%% Read metafile
metafile = 'wtk_site_metadata.csv';
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
data = struct(...
    'nI', nI,...
    'nT', nT,...
    'sid', sid,...
    'cap', cap,...
    'yyyy', yyyy,...
    'mm', mm,...
    'dd', dd,...
    'hh', hh,...
    'mi', mi,...
    'va', va,...
    'xa', xa,...
    'xf', xf,...
    'xa_MW', xa_MW,...
    'xf_MW', xf_MW...
    );
end