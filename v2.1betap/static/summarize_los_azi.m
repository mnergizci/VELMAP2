function track_summary = summarize_los_azi(data)
% summarize_los_azi - Compute mean LOS and azimuth for a stack (insar or sboi)
% 
% Usage:
%   summarize_los_azi(insar)
%   summarize_los_azi(sboi)
%
% Input:
%   data - struct array with fields: proc.ID, los, azi
%
% Output:
%   track_summary - struct array with ID, mean_los, mean_azi

numTracks = length(data);
track_summary = struct('ID', {}, 'mean_los', {}, 'mean_azi', {});

for i = 1:numTracks
    try
        ID = data(i).proc.ID;
        los = data(i).los;
        azi = data(i).azi;

        mean_los = mean(los(:), 'omitnan');
        mean_azi = mean(azi(:), 'omitnan');

        track_summary(end+1) = struct( ...
            'ID', ID, ...
            'mean_los', mean_los, ...
            'mean_azi', mean_azi);
    catch ME
        warning('Skipping entry %d due to missing fields: %s', i, ME.message);
    end
end

% Display results
fprintf('Track ID\t\tMean LOS\tMean Azimuth\n');
for i = 1:length(track_summary)
    fprintf('%s\t%.4f\t\t%.4f\n', ...
        track_summary(i).ID, ...
        track_summary(i).mean_los, ...
        track_summary(i).mean_azi);
end

end
