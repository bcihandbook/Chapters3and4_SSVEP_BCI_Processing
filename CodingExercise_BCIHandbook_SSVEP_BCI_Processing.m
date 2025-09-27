clear all;close all;clc;

%% Add bdfmatlab for loading data
addpath(genpath('.\bdfmatlab')); % Change to your own path

%% Useful init
load('chanlocs64.mat'); % Biosemi 64-channel map
SamplingRate = 2048; % Hz

%% Butterworth filter
[b, a] = butter(4, [2 30]/(SamplingRate/2), 'bandpass');

%% Find all BDF data
DataPath = ''; % Fill in the path to the data location on your machine
BDFfiles = dir([DataPath '/*.bdf']);
BDFfiles = {BDFfiles.name};


%% Load data, extract trials and labels
Trials = [];
Labels = [];
trial = 0;
for f=1:length(BDFfiles)

    header = readbdfheader([DataPath '/' BDFfiles{f}]);
    data = readbdfdata(header);
    data = data';
    
    %% Keep trigger channel
    trigger = data(:,end);
    
    %% Keep EEG channels
    ChannelLabels = {header.Channel.Label};
    ChannelLabels = ChannelLabels(1:64);
    data = data(:,[1:64]);
     
    %% Get events
    POS = gettrigger(trigger, 305);
    TYP = trigger(gettrigger(trigger,305));
    
    %% Trial extraction
    for event=1:length(TYP)
        if(TYP(event) ~= 768)
            trial = trial + 1;
            Trials(trial,:,:) = data(POS(event):POS(event)+SamplingRate*5-1,:);
            Labels(trial) = TYP(event);
        end
    end
    
end

montage = channels2montage(ChannelLabels);
lapplacianMatrix = montage2laplacian(montage);

%% Pre-processing and feature extraction for each trial
for tr=1:size(Trials,1)

    %% Band-pass filtering
    processedtrial = filtfilt(b,a,squeeze(Trials(tr,:,:)));
   
    %% DC removal (is it really needed after band-pass filtering?)
    processedtrial = removeDC(processedtrial);
    
    %% Spatial fitlering
    % CAR
    %processedtrial = car(processedtrial);
    
    %Laplacian
    processedtrial = processedtrial*lapplacianMatrix;
    
    % Keep a copy of the pre-processed trials to be used later on with the
    % CCA-based decoding method 
    pTrials(tr,:,:) = processedtrial;
    
    %% Extract spectral features for each channel
    for ch=1:size(processedtrial,2)
        ch
        % PSD
        [psd(tr,ch,:) psdfreqs] = pwelch(processedtrial(:,ch), SamplingRate*2,...
            SamplingRate*0.25 , [], SamplingRate);
        
        %% FFT (note: spectrogram will do FFT over time)
        [spectrumfft(tr,ch,:),~, fftfreqs] = bpfft(processedtrial(:,ch), SamplingRate);
        
        % Continuous Wavelet Transform (CWT) 
        % Note: if you want only the spectrum and not a time-frequency
        % decomposition, this kind of defeats the purpose of using wavelets, spectrogram, etc.)
        %[cfs, waveletfreq] = cwt(processedtrial(:,ch),'morse',SamplingRate);
        %spectrumcwt(tr,ch,:) = mean(abs(cfs),2);
        
        % There are million other ways to get a spectrum estimate
        % - Filter bank + squaring of the signal
        % - Multitapers
        % - And many others
    end
end
load('WaveletResults.mat'); % Load pre-saved CWT cause it takes ages

% Per-class spectra
ChannelName = 'O1'; % Change to view other channels
ChannelIndex = find(strcmp(ChannelLabels,ChannelName));
UniqueLabels = unique(Labels);
for c=1:length(UniqueLabels)
    GApsd(c,:) = squeeze(mean(psd(Labels==UniqueLabels(c),ChannelIndex,:)));
    GAfft(c,:) = squeeze(mean(spectrumfft(Labels==UniqueLabels(c),ChannelIndex,:)));
    GAcwt(c,:) = squeeze(mean(spectrumcwt(Labels==UniqueLabels(c),ChannelIndex,:)));
    
    STDpsd(c,:) = squeeze(std(psd(Labels==UniqueLabels(c),ChannelIndex,:)));
    STDfft(c,:) = squeeze(std(spectrumfft(Labels==UniqueLabels(c),ChannelIndex,:)));
    STDcwt(c,:) = squeeze(std(spectrumcwt(Labels==UniqueLabels(c),ChannelIndex,:)));
end

% PSD
figure();
lof = find(psdfreqs==1);
hif = find(psdfreqs==40);
plot(psdfreqs(lof:hif),GApsd(:,lof:hif),'LineWidth',2);
xlabel('Frequency [Hz]','FontSize',20);
ylabel('PSD Bandpower','FontSize',20);
title(['Channel ' ChannelName],'FontSize',20)
legend({'7.5 Hz','9 Hz','10 Hz','12 Hz','15 Hz','20 Hz'});

% FFT
figure();
lof = find(fftfreqs==1);
hif = find(fftfreqs==40);
plot(fftfreqs(lof:hif),GAfft(:,lof:hif),'LineWidth',2);
xlabel('Frequency [Hz]','FontSize',20);
ylabel('FFT Bandpower','FontSize',20);
title(['Channel ' ChannelName],'FontSize',20)
legend({'7.5 Hz','9 Hz','10 Hz','12 Hz','15 Hz','20 Hz'});

% CWT
figure();
%lof = find(waveletfreq >= 1);lof = lof(end);
%hif = find(waveletfreq<=40);
plot(waveletfreq(end-50:end),GAcwt(:,end-50:end),'LineWidth',2);
xlabel('Frequency [Hz]','FontSize',20);
ylabel('FFT Bandpower','FontSize',20);
title(['Channel ' ChannelName],'FontSize',20)
legend({'7.5 Hz','9 Hz','10 Hz','12 Hz','15 Hz','20 Hz'});

%% GAs with std info (example with PSD and FFT)
Colors = {'k','r','b','m','c','g'};

% PSD
figure();
lof = find(psdfreqs==1);
hif = find(psdfreqs==40);
% Add std info
for c=1:length(UniqueLabels)
    h(c) = plot(psdfreqs(lof:hif),GApsd(c,lof:hif),Colors{c},'LineWidth',2);hold on;
    shadedErrorBar(psdfreqs(lof:hif),GApsd(c,lof:hif),STDpsd(c,lof:hif),Colors{c},1);
end
hold off;
xlabel('Frequency [Hz]','FontSize',20);
ylabel('PSD Bandpower','FontSize',20);
title(['Channel ' ChannelName],'FontSize',20)
legend(h,{'7.5 Hz','9 Hz','10 Hz','12 Hz','15 Hz','20 Hz'});

% FFT
figure();
lof = find(fftfreqs==1);
hif = find(fftfreqs==40);
% Add std info
for c=1:length(UniqueLabels)
    h(c) = plot(fftfreqs(lof:hif),GAfft(c,lof:hif),Colors{c},'LineWidth',2);hold on;
    shadedErrorBar(fftfreqs(lof:hif),GAfft(c,lof:hif),STDfft(c,lof:hif),Colors{c},1);
end
hold off;
xlabel('Frequency [Hz]','FontSize',20);
ylabel('FFT Bandpower','FontSize',20);
title(['Channel ' ChannelName],'FontSize',20)
legend(h,{'7.5 Hz','9 Hz','10 Hz','12 Hz','15 Hz','20 Hz'});

%% Downsample as needed
TargetFreqs = [5:0.5:25];
IndPSD = ismember(psdfreqs, TargetFreqs);
IndFFT = ismember(fftfreqs, TargetFreqs);
psd = psd(:,:,IndPSD);
spectrumfft = spectrumfft(:,:,IndFFT);
psdfreqs = psdfreqs(IndPSD);
fftfreqs = fftfreqs(IndFFT);


%% Feature selection
% r^2, CVA or Fisher Score? Why?
for ch=1:length(ChannelLabels)
    for fr=1:length(TargetFreqs)
        % R2 
        r2psd(ch,fr) = rsquared(psd(:,ch,fr),Labels');
        r2fft(ch,fr) = rsquared(spectrumfft(:,ch,fr),Labels');
        
        % Fisher Score
        for c=1:length(UniqueLabels)
            tmplbl = nan(length(Labels),1);
            tmplbl(Labels==UniqueLabels(c))=1;
            tmplbl(Labels~=UniqueLabels(c))=2;
            fspsd(ch,fr,c) = fisherscore(psd(:,ch,fr),tmplbl);
            fsfft(ch,fr,c) = fisherscore(spectrumfft(:,ch,fr),tmplbl);
        end 
    end
end
afspsd = mean(fspsd,3);
afsfft = mean(fsfft,3);


figure();
subplot(2,2,1);imagesc(r2fft);
set(gca,'XTick',[1:1:length(TargetFreqs)]);
set(gca,'XTickLabel',TargetFreqs);
set(gca,'YTick',[1:1:length(ChannelLabels)]);
set(gca,'YTickLabel',ChannelLabels);
subplot(2,2,2);imagesc(r2psd);
set(gca,'XTick',[1:1:length(TargetFreqs)]);
set(gca,'XTickLabel',TargetFreqs);
set(gca,'YTick',[1:1:length(ChannelLabels)]);
set(gca,'YTickLabel',ChannelLabels);
subplot(2,2,3);imagesc(afsfft);
set(gca,'XTick',[1:1:length(TargetFreqs)]);
set(gca,'XTickLabel',TargetFreqs);
set(gca,'YTick',[1:1:length(ChannelLabels)]);
set(gca,'YTickLabel',ChannelLabels);
subplot(2,2,4);imagesc(afspsd);
set(gca,'XTick',[1:1:length(TargetFreqs)]);
set(gca,'XTickLabel',TargetFreqs);
set(gca,'YTick',[1:1:length(ChannelLabels)]);
set(gca,'YTickLabel',ChannelLabels);


%% Topoplot of max
figure();
subplot(1,2,1);
topoplot(max(afspsd'),chanlocs,'electrodes','labels','maplimits',[0 1]);colorbar;
subplot(1,2,2);
topoplot(max(afsfft'),chanlocs,'electrodes','labels','maplimits',[0 1]);colorbar;

%% Select type of features to use
Data = spectrumfft; % Change variable to use other types of features

%% Supervised classification
NFeatures = 8; % You can play around with different numbers of features
TrainingPrct = 0.5; % You can play around with different train/test split ratios
TrainingData = Data(1:round(TrainingPrct*size(Data,1)),:,:);
TrainingLabels = Labels(1:round(TrainingPrct*size(Data,1)));
TestingData = Data(round(TrainingPrct*size(Data,1))+1:end,:,:);
TestingLabels = Labels(round(TrainingPrct*size(Data,1))+1:end);
% Do feature selection only on training set
for ch=1:length(ChannelLabels)
    for fr=1:length(TargetFreqs)
        
        % Fisher Score
        for c=1:length(UniqueLabels)
            tmplbl = nan(length(TrainingLabels),1);
            tmplbl(TrainingLabels==UniqueLabels(c))=1;
            tmplbl(TrainingLabels~=UniqueLabels(c))=2;
            trainfs(ch,fr,c) = fisherscore(TrainingData(:,ch,fr),tmplbl);
        end 
    end
end
trainfs = mean(trainfs,3);

% Select N uncorrelated features
[~, BestFeatInd] = sort(trainfs(:),'descend');
[BestFeatChannel, BestFeatBand] = ind2sub(size(trainfs),BestFeatInd);

SelectedTrainingData = TrainingData(:, BestFeatChannel(1),BestFeatBand(1));
SelectedTestingData = TestingData(:, BestFeatChannel(1),BestFeatBand(1));
NSelected = 1;
sf = 2;
SelectedFeatChannel = [];
SelectedFeatBand = [];
while((NSelected < NFeatures) && (sf <= length(BestFeatInd)))
    
    candidateFeature = TrainingData(:, BestFeatChannel(sf),BestFeatBand(sf));
    corrprev = [];
    for prev=1:size(SelectedTrainingData,2)
        corrprev(prev) = corr(candidateFeature, SelectedTrainingData(:,prev));
    end
    if(max(corrprev) < 0.75) % You can change the threshold of acceptable correlation
        % Add this feature
        NSelected = NSelected + 1;
        SelectedTrainingData = [SelectedTrainingData candidateFeature];
        SelectedTestingData = [SelectedTestingData TestingData(:, BestFeatChannel(sf),BestFeatBand(sf))];
        SelectedFeatChannel = [SelectedFeatChannel ; BestFeatChannel(sf)];
        SelectedFeatBand = [SelectedFeatBand ; BestFeatBand(sf)];
    end
    
    % Try next feature
    sf = sf+1;
end

% Train and test with qda
QDAModel = fitcdiscr(SelectedTrainingData, TrainingLabels,'DiscrimType','quadratic'); % You can try linear (LDA) and diaglinear (Naive Bayes)
PredictedClass = predict(QDAModel, SelectedTestingData);
Accuracy = 100*sum(PredictedClass' == TestingLabels)/length(TestingLabels);
for c=1:length(UniqueLabels)
    ClassAccuracy(c) = 100*sum(PredictedClass(TestingLabels==UniqueLabels(c))' == UniqueLabels(c))/sum(TestingLabels==UniqueLabels(c));
end

% CCA -- No need for training/calibration (non-parametric model)
StimuliFrequencies = [7.5 9 10 12 15 20];
IdealSignals = {};
t = [0:1:size(pTrials,2)-1]'/SamplingRate; % Time vector
for stim=1:length(StimuliFrequencies)
    % For each stimulus, we add one ideal sine and one cosine, so as to be
    % able to caputre a 90 degree phase difference. More phases can be
    % added for the same reason. This can be extended for the first few
    % harmonics of each stimulus base frequency. Here, stimuli 7.5 and 10
    % Hz have a second harmonic that coincides with other stimuli, so,
    % adding the harmonics would not make sense.
    IdealSignals{stim} = [sin(2*pi*StimuliFrequencies(stim)*t) cos(2*pi*StimuliFrequencies(stim)*t)];
end

% Use a single channel for simplicity, without loss of generality (i.e. the
% method works similarly with many channels, but the prediction will fail
% if the additional channels carry no information)
ChannelName = 'O1';
ChannelIndex = find(strcmp(ChannelLabels,ChannelName));
for tr=1:size(pTrials,1)
    rho = [];
    for stim=1:length(StimuliFrequencies)
        [~,~, ccacor] = canoncorr(IdealSignals{stim}, squeeze(pTrials(tr,:,ChannelIndex)'));
        rho(stim) = max(ccacor);
    end
    [~, maxind] = max(rho);
    PredictedLabelsCCA(tr) = UniqueLabels(maxind);
end

% Note that here, we did not "learn" anything for the decoder, we did not
% estimate any parameters. Hence, this method is not in danger of
% overfitting, and we can thus use the whole dataset to test it, no need
% for cross-validation or training/testing split
AccuracyCCA = 100*sum(PredictedLabelsCCA == Labels)/length(Labels);
for c=1:length(UniqueLabels)
    ClassAccuracyCCA(c) = 100*sum(PredictedLabelsCCA(Labels==UniqueLabels(c))' == UniqueLabels(c))/sum(Labels==UniqueLabels(c));
end
% Conclusion: Using CCA (i.e., a method that needs no supervised learning
% and therefore no calibration at all, allowing direct closed-loop SSVEP 
% BCI without any boring, time-consuming training!) we get higher and more 
% balanced (across classes) classification accuracy than with a QDA machine
% learning classifier trained with supervised learning. Note that here we
% used a single channel and did not optimize the CCA method (could work
% better with downsmapling, using more channels, harmonics, etc.)