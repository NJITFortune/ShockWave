% batchEigenfreqs.m
% Batch processes Eigenmannia frequency tracking across sequential .mat files.
%
% For the first file, the user clicks to identify EOD frequencies for each
% channel.  For subsequent files, the mean frequencies from the previous file
% are passed as prefreqs so no clicking is required.
%
% Each .mat file is expected to contain:
%   data         - recorded signals [nChannels x nSamples]
%   sampleRate   - sample rate in Hz
%   channelNames - cell array of human-readable channel name strings
%   startTime    - datetime string for the start of the recording
%
% Output structure (results):
%   results(iFile).filename               - source .mat filename
%   results(iFile).startTime              - startTime from the file
%   results(iFile).channels(iChan).name  - channelNames entry for this channel
%   results(iFile).channels(iChan).pf    - tracked frequencies [nFish x nWindows]
%   results(iFile).channels(iChan).wtims - time stamps for each window

%% Select directory
dirPath = uigetdir(pwd, 'Select directory containing .mat data files');
if isequal(dirPath, 0)
    error('No directory selected. Exiting.');
end

%% Get sorted list of .mat files
fileList = dir(fullfile(dirPath, '*.mat'));
if isempty(fileList)
    error('No .mat files found in: %s', dirPath);
end

[~, sortIdx] = sort({fileList.name});
fileList = fileList(sortIdx);
nFiles = length(fileList);
fprintf('Found %d .mat file(s) in %s\n\n', nFiles, dirPath);

%% Process files
r = struct('filename', cell(nFiles,1), 'startTime', cell(nFiles,1), 'channels', cell(nFiles,1));
prefreqs = {};   % prefreqs{iChan} holds the mean freqs from the previous file
nTanks = [];  % determined from the first file

for iFile = 1:nFiles

    fname = fullfile(dirPath, fileList(iFile).name);
    fprintf('--- File %d/%d: %s ---\n', iFile, nFiles, fileList(iFile).name);

    % Load data
    loaded       = load(fname);
    data         = loaded.data;
    Fs           = loaded.sampleRate;
    tankNames    = loaded.channelNames;
    startTime    = loaded.startTime;

    % Determine channel count from the first file
    if iFile == 1
        nTanks = length(data(1,:));
        fprintf('Detected %d Tank(s).\n', nTanks);
        prefreqs = cell(nTanks, 1);   % all empty => user clicks on file 1
    end

    % Storage for this file's results
    chanResults(nTanks) = struct('name', [], 'pf', [], 'pa', [], 'wtims', []);

    for iChan = 1:nTanks
        chanName = tankNames{iChan};
        fprintf('  Tank %d/%d: %s ...\n', iChan, nTanks, chanName);

        chanData = data(:, iChan);

        [pf, pa, wtims] = findEigenfreqs(chanData, Fs, prefreqs{iChan});

        chanResults(iChan).name  = chanName;
        chanResults(iChan).pf    = pf;
        chanResults(iChan).pa    = pa;
        chanResults(iChan).wtims = wtims;

        % Mean frequency of each fish across the file -> prefreqs for next file
        prefreqs{iChan} = mean(pf, 2);
        
    end

    r(iFile).filename  = fileList(iFile).name;
    r(iFile).startTime = startTime;
    r(iFile).tankID  = chanResults;

    fprintf('  Done with %s.\n\n', fileList(iFile).name);
end

fprintf('Batch complete.  Processed %d file(s), %d channel(s) each.\n', nFiles, nTanks);
