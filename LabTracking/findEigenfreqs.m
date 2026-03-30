function [pf, pa, wtims] = findEigenfreqs(data, Fs, prefreqs)
% Reads a single channel of Eigenmannia data recorded in a single tank
% (data) and traces the EOD frequencies in the sample. pf is the list of
% EOD frequencies for each fish.  wtims is the time stamps for the middle
% of each FFT window used for estimating the EOD frequencies.
% Sadly, this script requires too much user clicking - to get the initial
% EOD frequencies and to correct errors.
% Frequencies are reported at the middle of the time window over which the
% long FFT was taken.
% This function has two embedded functions, fftMaker and getpeaks.
% fftMaker calculates the power spectrum of a sample.
% getpeaks uses fftMaker to choose 'new' peak frequencies based on previous peaks.

%% User changeable settings
freqRange = [250 550];           % Frequency range for Eigenmannia
nFFT = 32*1024;                  % Large FFT window for precise freq estimation
stepsize = ceil(nFFT * 0.05);    % value of 0.05 is 95% overlap

%% Pre calculations
% Make a time series for the data
tim = 1/Fs:1/Fs:length(data)/Fs;

% These are the centers of the time windows for which the frequencies of
% the fish will be calculated. wtims are the actual times, widxs are the
% indicies from tim.
wtims = tim(1+(nFFT/2):stepsize:(length(data)-(nFFT/2))-1);
% This loop is slow - there has to be a better way!  
%    for j=length(wtims):-1:1
%        widxs(j) = find(tim == wtims(j));
%    end
    widxs = 1+(nFFT/2) : stepsize : (length(data)-(nFFT/2))-1;                                                    

%% Interface with the user to get the initial frequencies

% Get clicks for each Eigen frequency

figure(1); clf; specgram(data, nFFT/2, Fs, [], floor(0.80*(nFFT/2))); ylim(freqRange);

if isempty(prefreqs)

    [startTim, userFreqs] = ginput();
    
    startTim = mean(startTim);
    
    % The user may click at either end or in the middle. This handles that.
    if startTim < wtims(1)
        startTim = wtims(1);
        direction = 1;
    elseif startTim > wtims(end)
        startTim = wtims(end);
        direction = 2; 
    else
        startTim = wtims(find(wtims >= startTim, 1));
        direction = 3;
    end
    
    startWidx = find(wtims == startTim);

else
    userFreqs = prefreqs;
    startWidx = 1;
    direction = 1;

    for j = 1:length(prefreqs)
        figure(1); hold on; plot([0, 0.5], [prefreqs(j), prefreqs(j)]);
    end

    if length(prefreqs) > 1
        % Ask if we should reclick

    end
end

pf(1:length(userFreqs),length(wtims)) = zeros(1,length(userFreqs));
pa(1:length(userFreqs),length(wtims)) = zeros(1,length(userFreqs));

newFreqs1 = userFreqs; % Lazy coding - newFreqs1 gets updated for forward analysis.
newFreqs2 = userFreqs; % Lazy stupid coding - newFreqs2 gets update for the backward analysis. Embarassing.

%% Automated frequency tracking

% Go forward young computer - forwards in time from where the user specified the EOD frequencies
if direction ~= 2
    for j = startWidx:length(wtims)

       curWindowIDX = widxs(j)-(nFFT/2):widxs(j)+(nFFT/2);    % This is the indicies for the window of data to analyze
       [pf(:,j), pa(:,j)] = getpeaks(data(curWindowIDX), Fs, newFreqs1); % Get the current peaks based on previous peaks
       newFreqs1 = pf(:,j);                                   % Assign the values to the our data out.

    end
end

% Go backward young computer - now go backwards in time from where the user specified the EOD frequencies.
if direction ~= 1
    for j = startWidx:-1:1
        
       curWindowIDX = widxs(j)-(nFFT/2):widxs(j)+(nFFT/2);
       [pf(:,j), pa(:,j)] = getpeaks(data(curWindowIDX), Fs, newFreqs2);
       newFreqs2 = pf(:,j);

    end
end

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
    while keepGoingWhileOne == 1 && curIDX(end) <= length(wtims)
        if length(unique(pf(:,curIDX(end)))) == length(userFreqs)
            keepGoingWhileOne = 2;
        end
        curIDX(end+1) = curIDX(end) + 1;
    end
    curIDX = curIDX(curIDX <= length(wtims));

    if isempty(prefreqs)
        f = fftMaker(data(widxs(1)-(nFFT/2):widxs(1)+(nFFT/2)), Fs, 3);
        ff.freqs = f.fftfreq(f.fftfreq > 250 & f.fftfreq < 650);
        ff.data = f.fftdata(f.fftfreq > 250 & f.fftfreq < 650);

        figure(2); clf; plot(ff.freqs, ff.data, 'k');
        hold on; text(freqRange(1)+40,0, [num2str(length(userFreqs)), ' Start'], 'FontSize', 24);

        [userNewFreqs, ~] = ginput(length(userFreqs));
    else
        fprintf('    Fix #1 (start): using prefreqs as seed (no clicks needed).\n');
        userNewFreqs = userFreqs;
    end

    for j=curIDX
        curWindowIDX = widxs(j)-(nFFT/2):widxs(j)+(nFFT/2);
        [pf(:,j), pa(:,j)]= getpeaks(data(curWindowIDX), Fs, userNewFreqs);
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

    while keepGoingWhileEndLen == length(wtims) && rucIDX(end) >= 1
        if length(unique(pf(:,rucIDX(end)))) == length(userFreqs)
            keepGoingWhileEndLen = 2;
        end
        rucIDX(end+1) = rucIDX(end) - 1;
    end
    rucIDX = rucIDX(rucIDX >= 1);

    if isempty(prefreqs)
        f = fftMaker(data(widxs(end)-(nFFT/2):widxs(end)+(nFFT/2)), Fs, 3);
        ff.freqs = f.fftfreq(f.fftfreq > 250 & f.fftfreq < 650);
        ff.data = f.fftdata(f.fftfreq > 250 & f.fftfreq < 650);

        figure(2); clf; plot(ff.freqs, ff.data, 'k');
        hold on; text(freqRange(1)+20,0, [num2str(length(userFreqs)), ' End'], 'FontSize', 24);

        [userNewFreqs, ~] = ginput(length(userFreqs));
    else
        fprintf('    Fix #2 (end): using prefreqs as seed (no clicks needed).\n');
        userNewFreqs = userFreqs;
    end

    for j=rucIDX
           curWindowIDX = widxs(j)-(nFFT/2):widxs(j)+(nFFT/2);
           [pf(:,j), pa(:,j)] = getpeaks(data(curWindowIDX), Fs, userNewFreqs);
           userNewFreqs = pf(:,j);
    end

end

%% Plot the result

figure(1); clf; specgram(data, nFFT, Fs, [], floor(0.80*nFFT)); ylim(freqRange);
colormap('HOT'); 
hold on;
for j=1:length(userFreqs)

    plot(wtims, pf(j,:), '.', 'MarkerSize', 12);

end

ylim([min(pf(1,:))-20, max(pf(end,:))+20]);

% Fix #3 - The two previous fixes don't solve crossovers that occur in the middle
% of the recording.  Here the user clicks the location, then the EOD frequencies,
% and the code retraces.  Then we have to align the indices. 

end

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
 
function [peakfreqs, peakamps] = getpeaks(snip, samplerate, prefreqs)

    m = fftMaker(snip, samplerate, 3);

    for j=length(prefreqs):-1:1

        % Fundamental: search ±1 Hz around previous estimate
        f1idx = m.fftfreq > prefreqs(j)-1 & m.fftfreq < prefreqs(j)+1;
        [a1, i1] = max(m.fftdata(f1idx));
        f1freqs = m.fftfreq(f1idx);
        f1 = f1freqs(i1);

        % 1st harmonic: search ±2 Hz around 2x previous estimate (window
        % doubled because absolute drift is twice as large at the harmonic)
        h1idx = m.fftfreq > 2*prefreqs(j)-2 & m.fftfreq < 2*prefreqs(j)+2;
        if any(h1idx)
            [a2, i2] = max(m.fftdata(h1idx));
            h1freqs = m.fftfreq(h1idx);
            f2 = h1freqs(i2) / 2;   % convert harmonic peak back to fundamental
            peakfreqs(j) = (f1*a1 + f2*a2) / (a1 + a2);  % amplitude-weighted mean
        else
            peakfreqs(j) = f1;       % harmonic out of range, use fundamental only
        end

        peakamps(j) = a1;            % report fundamental amplitude

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
