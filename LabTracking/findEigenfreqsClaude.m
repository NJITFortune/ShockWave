function [pf, wtims] = findEigenfreqs(data, Fs)
% Reads a single channel of Eigenmannia data recorded in a single tank
% (data) and traces the EOD frequencies in the sample. pf is the list of
% EOD frequencies for each fish.  wtims is the time stamps for the middle
% of each FFT window used for estimating the EOD frequencies.
% Initial EOD frequencies and the starting window are detected automatically.
% Frequencies are reported at the middle of the time window over which the
% long FFT was taken.
% Embedded functions: fftMaker, getpeaks, removeHarmonics.
% fftMaker computes the power spectrum of a snippet.
% getpeaks finds updated peak frequencies near a set of prior frequencies.
% removeHarmonics discards peaks that are integer multiples of a lower peak.

%% User changeable settings
freqRange = [250 650];           % Frequency range for Eigenmannia
nFFT = 16*1024;                  % Large FFT window for precise freq estimation
stepsize = ceil(nFFT * 0.05);    % value of 0.05 is 95% overlap

%% Pre calculations
% Make a time series for the data
tim = 1/Fs:1/Fs:length(data)/Fs;

% These are the centers of the time windows for which the frequencies of
% the fish will be calculated. wtims are the actual times, widxs are the
% indicies from tim.
wtims = tim(1+(nFFT/2):stepsize:(length(data)-(nFFT/2))-1);
widxs = 1+(nFFT/2) : stepsize : (length(data)-(nFFT/2))-1;

%% Automatically detect fish count and find the best starting window

% Get frequency axis (identical for every window since all windows are nFFT+1 samples)
f0 = fftMaker(data(widxs(1)-(nFFT/2):widxs(1)+(nFFT/2)), Fs, 3);
inRange    = f0.fftfreq >= freqRange(1) & f0.fftfreq <= freqRange(2);
rangeFreqs = f0.fftfreq(inRange)';                     % column vector
freqStep   = rangeFreqs(2) - rangeFreqs(1);
minSepSamp = round(3 / freqStep);                      % 3 Hz minimum separation

% --- Step 1: Estimate fish count from a median spectrum ---
% Averaging over many windows suppresses transient events; the median is
% more robust to the occasional frequency crossing or noise burst.
scanIdxs = unique(round(linspace(1, length(wtims), min(50, length(wtims)))));
specMat   = zeros(length(scanIdxs), length(rangeFreqs));
for ki = 1:length(scanIdxs)
    fk = fftMaker(data(widxs(scanIdxs(ki))-(nFFT/2):widxs(scanIdxs(ki))+(nFFT/2)), Fs, 3);
    specMat(ki,:) = fk.fftdata(inRange);
end
medSpec = median(specMat, 1);

[~, locs] = findpeaks(medSpec, ...
    'MinPeakProminence', 10 * median(medSpec), ...
    'MinPeakDistance',   minSepSamp);
candidateFreqs = rangeFreqs(locs);

% Remove harmonics: discard any peak within 3 Hz of an integer multiple
% (2x–5x) of a lower-frequency peak. Fish at ~300 Hz produce a 2nd harmonic
% at ~600 Hz that would otherwise be counted as a separate fish.
candidateFreqs = removeHarmonics(candidateFreqs, 3);
nFish = length(candidateFreqs);
fprintf('Estimated %d fish (median-spectrum peaks: %s Hz)\n', ...
    nFish, num2str(candidateFreqs', '%.1f '));

% --- Step 2: Find the window with the most-separated peaks ---
% Scan ~100 windows and keep the one where exactly nFish peaks are found
% and the minimum pairwise gap between them is largest. This naturally
% avoids frequency crossings, where two fish overlap and the gap collapses.
scanIdxs2 = unique(round(linspace(1, length(wtims), min(100, length(wtims)))));
bestSep   = 0;
bestWidx  = round(length(wtims) / 2);     % fallback to middle
for ki = 1:length(scanIdxs2)
    k  = scanIdxs2(ki);
    fk = fftMaker(data(widxs(k)-(nFFT/2):widxs(k)+(nFFT/2)), Fs, 3);
    rData = fk.fftdata(inRange);
    [~, locsK] = findpeaks(rData, ...
        'MinPeakProminence', 10 * median(rData), ...
        'MinPeakDistance',   minSepSamp);
    peaksK = removeHarmonics(rangeFreqs(locsK), 3);
    if length(peaksK) == nFish
        minSep = min(diff(sort(peaksK)));
        if minSep > bestSep
            bestSep  = minSep;
            bestWidx = k;
        end
    end
end

startWidx = bestWidx;
startTim  = wtims(startWidx);
direction = 3;    % always trace both directions from the chosen window

% Get the precise starting frequencies from the best window
fBest  = fftMaker(data(widxs(startWidx)-(nFFT/2):widxs(startWidx)+(nFFT/2)), Fs, 3);
rDataB = fBest.fftdata(inRange);
[~, locsB] = findpeaks(rDataB, ...
    'MinPeakProminence', 10 * median(rDataB), ...
    'MinPeakDistance',   minSepSamp);
userFreqs = removeHarmonics(rangeFreqs(locsB), 3);

fprintf('Starting at t = %.1f s  (min fish separation = %.1f Hz)\n', startTim, bestSep);

% Show spectrogram with detected frequencies and start time marked
figure(1); clf; specgram(data, nFFT/2, Fs, [], floor(0.80*(nFFT/2))); ylim(freqRange);
hold on;
for j = 1:length(userFreqs)
    yline(userFreqs(j), 'r--', sprintf('%.1f Hz', userFreqs(j)), ...
        'LabelVerticalAlignment', 'bottom');
end
xline(startTim, 'g-', sprintf('start (sep=%.1fHz)', bestSep), ...
    'LabelVerticalAlignment', 'top');
title(sprintf('Auto-detected %d fish, start t=%.1fs', nFish, startTim));
drawnow;

pf(1:length(userFreqs), length(wtims)) = zeros(1, length(userFreqs));

newFreqs1 = userFreqs;
newFreqs2 = userFreqs;

%% Automated frequency tracking

% Go forward young computer - forwards in time from where the user specified the EOD frequencies
if direction ~= 2
    for j = startWidx:length(wtims)
        
       curWindowIDX = widxs(j)-(nFFT/2):widxs(j)+(nFFT/2);    % This is the indicies for the window of data to analyze
       pf(:,j) = getpeaks(data(curWindowIDX), Fs, newFreqs1); % Get the current peaks based on previous peaks
       newFreqs1 = pf(:,j);                                   % Assign the values to the our data out.

    end
end

% Go backward young computer - now go backwards in time from where the user specified the EOD frequencies.
if direction ~= 1
    for j = startWidx:-1:1
        
       curWindowIDX = widxs(j)-(nFFT/2):widxs(j)+(nFFT/2);
       pf(:,j) = getpeaks(data(curWindowIDX), Fs, newFreqs2);
       newFreqs2 = pf(:,j);

    end
end

% Crossover fixes only apply when there are 2+ fish.
if nFish > 1

% If EOD frequencies cross, it is inevitable that the tracking of two
% "channels" will merge at the crossing point and stay together for the
% rest of the sample.

% Fix #1 - see if our first sample has fewer unique EODs than expected.
% What we do is see if there are too few EOD frequencies at the beginning 
% of the recording, and then continue forwards until we get to a time sample 
% with the correct number.

keepGoingWhileOne = 1; % Again, I am terrible at coding.

if length(unique(pf(:,1))) < length(userFreqs)

    curIDX(1) = keepGoingWhileOne;
    while keepGoingWhileOne == 1
        if length(unique(pf(:,curIDX(end)))) == length(userFreqs)
            keepGoingWhileOne = 2;
        end
        curIDX(end+1) = curIDX(end) + 1;
    end

    f = fftMaker(data(widxs(1)-(nFFT/2):widxs(1)+(nFFT/2)), Fs, 3);
    ff.freqs = f.fftfreq(f.fftfreq > 250 & f.fftfreq < 650);
    ff.data = f.fftdata(f.fftfreq > 250 & f.fftfreq < 650);

    figure(2); clf; plot(ff.freqs, ff.data, 'k'); 

    [userNewFreqs, ~] = ginput();

    for j=curIDX
        curWindowIDX = widxs(j)-(nFFT/2):widxs(j)+(nFFT/2);
        pf(:,j) = getpeaks(data(curWindowIDX), Fs, userNewFreqs);
        userNewFreqs = pf(:,j);
    end

end

% Fix #2 - see if our last sample has fewer unique EODs than expected.
% What we do is see if there are too few EDO frequencies at the END 
% of the recording, and then continue backwards until we get to a time sample 
% with the correct number.

keepGoingWhileEndLen = length(wtims); % This is embarassing, isn't it?

if length(unique(pf(:,end))) < length(userFreqs)

    rucIDX(1) = keepGoingWhileEndLen;

    while keepGoingWhileEndLen == length(wtims)
        if length(unique(pf(:,rucIDX(end)))) == length(userFreqs)
            keepGoingWhileEndLen = 2;
        end
        rucIDX(end+1) = rucIDX(end) - 1;
    end

    f = fftMaker(data(widxs(end)-(nFFT/2):widxs(end)+(nFFT/2)), Fs, 3);
    ff.freqs = f.fftfreq(f.fftfreq > 250 & f.fftfreq < 650);
    ff.data = f.fftdata(f.fftfreq > 250 & f.fftfreq < 650);

    figure(2); clf; plot(ff.freqs, ff.data, 'k'); 

    [userNewFreqs, ~] = ginput();

    for j=rucIDX
           curWindowIDX = widxs(j)-(nFFT/2):widxs(j)+(nFFT/2);
           pf(:,j) = getpeaks(data(curWindowIDX), Fs, userNewFreqs);
           userNewFreqs = pf(:,j);
    end

end

end % nFish > 1

%% Plot the result

figure(1); clf; specgram(data, nFFT, Fs, [], floor(0.80*nFFT)); ylim(freqRange);
hold on;
for j=1:length(userFreqs)

    plot(wtims, pf(j,:), '.', 'MarkerSize', 16);

end

ylim([min(pf(1,:))-20, max(pf(end,:))+20]);

% Fix #3 - The two previous fixes don't solve crossovers that occur in the middle
% of the recording.  Here the user clicks the location, then the EOD frequencies,
% and the code retraces.  Then we have to align the indices. 




end

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%

function freqs = removeHarmonics(freqs, tol)
% Remove frequencies that are within tol Hz of an integer multiple (2x-5x)
% of a lower-frequency candidate. Input/output are column vectors.
    freqs = sort(freqs(:));
    keep  = true(size(freqs));
    for k = 1:length(freqs)
        for m = 1:k-1
            if keep(m)
                for n = 2:5
                    if abs(freqs(k) - n * freqs(m)) < tol
                        keep(k) = false;
                        break;
                    end
                end
            end
            if ~keep(k), break; end
        end
    end
    freqs = freqs(keep);
end

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%

function peakfreqs = getpeaks(snip, samplerate, prefreqs)

    m = fftMaker(snip, samplerate, 3);

    for j=length(prefreqs):-1:1
        freqfreqIDX = find(m.fftfreq > prefreqs(j)-1 & m.fftfreq < prefreqs(j)+1);
        [~,mIDX] = max(m.fftdata(freqfreqIDX));
        peakfreqs(j) = m.fftfreq(freqfreqIDX(mIDX));
    end

end

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 

function out = fftMaker(data, Fs, smoothwindow)
% Compute the FFT (Fast Fourier Transform)
% out = fftmachine(data, Fs, smoothwindow);
% Where out is a strucutre with fftfreq and fftdata
% The smoothwindow is for a medfilt1 low-pass filtering
% of the fft data itself.  This should generally be low and
% odd, 9 or less.

if nargin < 2
	fprintf('Usage: fftmachine(data, SampleRate, [window])');
	exit 0
end

L = length(data);

NFFT = 2^nextpow2(L); % Next power of 2 from length of the data

% NFFT = 1024*2;

fftdata = fft(data,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);

% fftdata = fft(data);

% We use only half of the data, hence fftdata(1:round(end/2));
% And we take the absolute value of the real component and filter
% that so that it is smooth

out.fftdata = 2*abs(fftdata(1:(NFFT/2)+1));
% out.fftdata = abs(real(fftdata(1:round(end/2))));

if nargin == 3
	out.fftdata = medfilt1( out.fftdata, smoothwindow);
end

% Now we need to generate the X values - which are the frequencies

%stepsize = Fs/round(length(data));
%out.fftfreq = stepsize:stepsize:Fs/2;
out.fftfreq = Fs/2*linspace(0,1,NFFT/2+1);

% Sometimes the rounding makes it so that the lengths of the
% data and the frequency values are off by one.  Let us correct that.

minlen = min([length(out.fftfreq) length(out.fftdata)]);
out.fftfreq = out.fftfreq(1:minlen);
out.fftdata = out.fftdata(1:minlen);

end