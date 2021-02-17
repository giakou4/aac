function AACSeq2 = AACoder2(fNameIn)
% function AACSeq2 = AACoder2(fNameIn)
% AACODER2 
%
% Inputs:
% fnameIN    : Name of file in .wav format. Sound must have 2 channels with
%              sampling frequency equals to 48kHz
%
% Outputs:
% AACSeq2    : Struct [K x 1] where K is the number of frames coded.
%              Contains:
%              * AACSeq2(i).frameType
%              * AACSeq2(i).winType
%              * AACSeq2(i).chl.TNScoeffs - quantized TNS coefficients of
%                                           left channel
%              * AACSeq2(i).chr.TNScoeffs - quantized TNS coefficients of
%                                           right channel
%              * AACSeq2(i).chl.frameF - MDCT coefficient of left channel
%                                      after TNS
%              * AACSeq2(i).chr.frameF - MDCT coefficient of right channel
%                                      after TNS
%
% Example:
% AACSeq2 = AACoder2('LicorDeCalandraca.wav');

% Add path of Level 1
warning off % I dont want to spamm warning for non invertible R in TNS
addpath([pwd '\..\Level 1'])

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
lineLength = fprintf('%2.0f / %4.0f',1,K-1);

% And so it begins 
for k=1:(K-1) % K-1 so we can get frameTnext
    if mod(k,10)==0
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf('%2.0f / %4.0f',k,K-1);
    end
    
    % Find current, next and 2 previous Frames in Time domain
    frameT = x((1024*k - 1023):(1024*k+1024),:); 
    frameTnext = x((1024*k + 1):(1024*k+2048),:); 
    if (1024*k-2047) < 1
        frameTprevious1 = zeros(1024,2);
    else
        frameTprevious1 = x((1024*k-2047):(1024*k),:);
    end
    if (1024*k-3071)< 1
        frameTprevious2 = zeros(1024,2);
    else
        frameTprevious2 = x((1024*k-3071):(1024*k-1024),:);
    end
    
    % Apply SSC to get frame type
    AACSeq2(k).frameType = SSC(frameT, frameTnext, frameTpreviousType);
    frameTpreviousType = AACSeq2(k).frameType;
    AACSeq2(k).winType = winType;
    
    % Apply Filterbank to get Frame in Frequency domain
    frameF = filterbank(frameT, AACSeq2(k).frameType, winType);
    
    % Apply TNS to get new coeffs in Frequency domain and add to AACSeq2
    if AACSeq2(k).frameType == "ESH"
        [AACSeq2(k).chl.frameF, AACSeq2(k).chl.TNScoeffs] = ...
                 TNS(reshape(frameF(:,1,:),[128 8]), AACSeq2(k).frameType);
        [AACSeq2(k).chr.frameF, AACSeq2(k).chr.TNScoeffs] = ...
                 TNS(reshape(frameF(:,2,:),[128 8]), AACSeq2(k).frameType);
            
    else
        [AACSeq2(k).chl.frameF, AACSeq2(k).chl.TNScoeffs] = ...
                                TNS(frameF(:,1), AACSeq2(k).frameType);
        [AACSeq2(k).chr.frameF, AACSeq2(k).chr.TNScoeffs] = ...
                                TNS(frameF(:,2), AACSeq2(k).frameType);
    end
    
end  

% Last printing and return
fprintf(repmat('\b',1,lineLength))
lineLength = fprintf('%2.0f / %4.0f',K-1,K-1);
fprintf('\n');

end

