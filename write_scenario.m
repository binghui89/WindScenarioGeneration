function write_scenario(sid, xa, xf, xs, header, dir)
% Write scenario data to csv file with given headers
% sid: nI by 1, xa and xf dim: nT by nI, xs dim: nT by nI by nS
% header: cell vector of (nS + 2) by 1, including xa and xf
if nargin == 5 % For compatibility with old code.
    dirwork = 'scenario_data';
else
    dirwork = dir;
end
dirhome = pwd;
commaHeader = [header;repmat({','},1,numel(header))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

nI = numel(sid);
nT = size(xa, 1);
nS = size(xs, 3);

if ~isdir(dirwork)
    mkdir(dirwork);
end
cd(dirwork);
for i = 1: nI
    %write header to file
    fname = strcat(int2str(sid(i)), '.csv');
    fid = fopen(fname,'w'); 
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
    
    %write data to end of file
    yourdata = nan(nT, nS+2); % Include actual and deterministic forecast
    yourdata = [xa(:, i), xf(:, i), squeeze(xs(:, i, :))];
    dlmwrite(fname,yourdata,'-append');
end
cd(dirhome);
fprintf('Scenario data written to .\\%s.\n', dirwork);

end