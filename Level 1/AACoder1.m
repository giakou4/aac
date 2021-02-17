function AACSeq1 = AACoder1(fNameIn)
% function AACSeq1 = AACoder1(fNameIn)
% AACODER1
%
% Inputs:
% fnameIN    : Name of file in .wav format. Sound must have 2 channels with
%              sampling frequency equals to 48kHz
%
% Outputs:
% AACSeq1    : Struct [K x 1] where K is the number of frames coded.
%              Contains:
%              * AACSeq1(i).frameType - OLS|LSS|LPPS|ESH as string
%              * AACSeq1(i).winType - SIN|KBD as string
%              * AACSeq1(i).chl.frameF - MDCT coefficient of left channel
%                                        [1024 x 1] for OLS|LSS|LPPS
%                                        [128 x 8] for ESH
%              * AACSeq1(i).chr.frameF - MDCT coefficient of right channel
%                                        [1024 x 1] for OLS|LSS|LPPS
%                                        [128 x 8] for ESH
%
% Example:
% AACSeq1 = AACoder1('LicorDeCalandraca.wav');

% Open audio file
[x,Fs] = audioread(fNameIn);
if Fs ~= 48000
    error('Sampling frequency is not 48.000Hz');
end

% Init and adjustments
N = length(x);
extraSamples = 1024 - mod(N, 1024); % for a full frame
x = [zeros(2048,2); x; zeros(2048+extraSamples,2)]; % zero padding 2048 samples and some extras
N = length(x);
K = floor(N/1024)-1; % find number of frames
winType = "SIN";
AACSeq1 = struct('frameType', num2cell(zeros(K-1,1)),...
                 'winType', num2cell(zeros(K-1,1)),...
                 'chl', struct('frameF',zeros(1024,1)),...
                 'chr', struct('frameF',zeros(1024,1)));
frameTpreviousType = "";


% Some printing    
fprintf('Iteration: '); 
lineLength = fprintf('%2.0f / %4.0f',1,(K-1));

% And so it begins  
for k=1:(K-1) % K-1 so we can get frameTnext => I throw last frame
    if mod(k,10)==0
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf('%2.0f / %4.0f',k,(K-1));
    end
    
    % Find current and next Frame in Time domain
    frameT = x((1024*k-1023):(1024*k+1024),:); 
    frameTnext = x((1024*k+1):(1024*k+2048),:); 
    
    % Apply SSC to get frame type
    AACSeq1(k).frameType = SSC(frameT, frameTnext, frameTpreviousType);
    frameTpreviousType = AACSeq1(k).frameType;
    AACSeq1(k).winType = winType;
    
    % Apply Filterbank to get Frame in Frequency domain and add to AACSeq1
    frameF = filterbank(frameT, AACSeq1(k).frameType, AACSeq1(k).winType);
    if AACSeq1(k).frameType == "ESH"
        AACSeq1(k).chl.frameF = reshape(frameF(:,1,:),[128 8]);
        AACSeq1(k).chr.frameF = reshape(frameF(:,2,:),[128 8]);
    else
        AACSeq1(k).chl.frameF = frameF(:,1);
        AACSeq1(k).chr.frameF = frameF(:,2);
    end
    
end

% Last printing and return
fprintf(repmat('\b',1,lineLength))
lineLength = fprintf('%2.0f / %4.0f',K-1,K-1);
fprintf('\n');
end

