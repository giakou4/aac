function AACSeq3 = AACoder3(fNameIn, fnameAACoded)
% function AACSeq3 = AACoder3(fNameIn)
% AACODER2 
% Inputs:
% fNameIn      : Name of file in .wav format. Sound must have 2 channels with
%                sampling frequency equals to 48kHz
% fnameAACoded : Name of .mat file to save AACSeq3
%
% Outputs:
% AACSeq3      : Struct [K x 1] where K is the number of frames coded.
%                Contains:
%                * AACSeq3(i).frameType
%                * AACSeq3(i).winType
%                * AACSeq3(i).chl.TNScoeffs - quantized TNS coefficients of
%                                           left channel
%                * AACSeq3(i).chr.TNScoeffs - quantized TNS coefficients of
%                                           right channel
%                * AACSeq3(i).chl.T - thresholds of psychoacoustic model
%                * AACSeq3(i).chr.T - thresholds of psychoacoustic model
%                * AACSeq3(i).chl.G - quantized global gains [8 x 1] for
%                                     ESH, [1 x 1] else
%                * AACSeq3(i).chr.G - quantized global gains [8 x 1] for 
%                                     ESH, [1 x 1] else
%                * AACSeq3(i).chl.sfc - Huffman coded sequence of sfc
%                                       ['101010']
%                * AACSeq3(i).chr.sfc - Huffman coded sequence of sfc
%                                       ['101010']
%                * AACSeq3(i).chl.stream - Huffman coded sequence of frame  
%                                          in frequency domain ['1010...']
%                * AACSeq3(i).chr.stream - Huffman coded sequence of frame
%                                          in frequency domain ['1010...']
%                * AACSeq3(i).chl.codebook - Huffman's codebook [int]
%                * AACSeq3(i).chr.codebook - Huffman's codebook [int]
%
% Example:
% AACSeq3 = AACoder3('LicorDeCalandraca.wav','LicorDeCalandraca_AACSeq3_demoAAC3.mat');

% Add path of Level 1 and 2
addpath([pwd '\..\Level 1'])
addpath([pwd '\..\Level 2'])

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
AACSeq3 = struct('frameType', num2cell(zeros(K-2,1)), ...
                 'winType', num2cell(zeros(K-2,1)), ...
                 'chl', struct('TNScoeffs', zeros(4,1), 'T', zeros(69,1), 'G', 0, 'sfc', '', 'stream', '', 'codebook', 0), ...
                 'chr', struct('TNScoeffs', zeros(4,1), 'T', zeros(69,1), 'G', 0, 'sfc', '', 'stream', '', 'codebook', 0));
frameTpreviousType = "";
huffLUT = loadLUT();

% Some printing   
fprintf('Iteration: '); 
lineLength = fprintf('%2.0f / %4.0f',1,K-1);

% And so it begins 
for k=1:(K-1) 
    if mod(k,10)==0
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf('%2.0f / %4.0f',k,K-1);
    end
    
    % Find current,previous, pre-previous and next Frame in Time domain
    frameT = x((1024*k - 1023):(1024*k+1024),:); 
    frameTnext = x((1024*k + 1):(1024*k+2048),:); 
    if (1024*k-2047) < 1
        frameTprevious1 = zeros(2048,2);
    else
        frameTprevious1 = x((1024*k-2047):(1024*k),:);
    end
    if (1024*k-3071)< 1
        frameTprevious2 = zeros(2048,2);
    else
        frameTprevious2 = x((1024*k-3071):(1024*k-1024),:);
    end
    
    % Apply SSC to get frame type
    AACSeq3(k).frameType = SSC(frameT, frameTnext, frameTpreviousType);
    frameTpreviousType = AACSeq3(k).frameType;
    AACSeq3(k).winType = winType;
    
    % Apply Filterbank to get Frame in Frequency domain
    frameF = filterbank(frameT, AACSeq3(k).frameType, winType); % [1024 x 2] or [128 x 2 x 8]
    
    if AACSeq3(k).frameType == "ESH"
         % Appy TNS to get transformed frameF
        [frameF(:,1,:), AACSeq3(k).chl.TNScoeffs] = ...
                              TNS(reshape(frameF(:,1,:),[128 8]), AACSeq3(k).frameType);
        [frameF(:,2,:), AACSeq3(k).chr.TNScoeffs] = ...
                              TNS(reshape(frameF(:,2,:),[128 8]), AACSeq3(k).frameType);        
        % Apply psychoacoustic model
        SMR = zeros(42,8,2);
        SMR(:,:,1) = psycho(frameT(:, 1), AACSeq3(k).frameType, frameTprevious1(:, 1), frameTprevious2(:, 1));
        SMR(:,:,2) = psycho(frameT(:, 2), AACSeq3(k).frameType, frameTprevious1(:, 2), frameTprevious2(:, 2));
        % Calculate T
        load('TableB219.mat');
        w_low = B219b(:,2);
        w_high = B219b(:,3);
        P = zeros(42,8);
        AACSeq3(k).chl.T = zeros(42,8); % reshape
        AACSeq3(k).chr.T = zeros(42,8); % reshape
        for i = 1:8
            for b = 1:42 
                AACSeq3(k).chl.T(b,i) = sum(frameF((w_low(b)+1):(w_high(b)+1),1,i).^2, 1);
                AACSeq3(k).chr.T(b,i) = sum(frameF((w_low(b)+1):(w_high(b)+1),2,i).^2, 1);
            end
        end
        AACSeq3(k).chl.T = AACSeq3(k).chl.T ./ SMR(:,:,1);
        AACSeq3(k).chr.T = AACSeq3(k).chr.T ./ SMR(:,:,2);
        % Apply Quantizer
        S = zeros(1024,2);
        sfc = zeros(42,8,2);
        AACSeq3(k).chl.G = zeros(8,1); % reshape
        AACSeq3(k).chr.G = zeros(8,1); % reshape
        [S(:,1), sfc(:,:,1), AACSeq3(k).chl.G] = AACquantizer(reshape(frameF(:,1,:),[128 8]), AACSeq3(k).frameType, SMR(:,:,1));
        [S(:,2), sfc(:,:,2), AACSeq3(k).chr.G] = AACquantizer(reshape(frameF(:,2,:),[128 8]), AACSeq3(k).frameType, SMR(:,:,2));
        % Apply Huffman
        [AACSeq3(k).chl.stream, AACSeq3(k).chl.codebook] = encodeHuff(S(:,1), huffLUT);
        [AACSeq3(k).chr.stream, AACSeq3(k).chr.codebook] = encodeHuff(S(:,2), huffLUT);
        [AACSeq3(k).chl.sfc, ~] = encodeHuff(reshape(sfc(:,:,1),[42*8 1]), huffLUT);
        [AACSeq3(k).chr.sfc, ~] = encodeHuff(reshape(sfc(:,:,2),[42*8 1]), huffLUT);   
    else
        % Appy TNS to get transformed frameF
        [frameF(:,1), AACSeq3(k).chl.TNScoeffs] = ...
                                TNS(frameF(:,1), AACSeq3(k).frameType);
        [frameF(:,2), AACSeq3(k).chr.TNScoeffs] = ...
                                TNS(frameF(:,2), AACSeq3(k).frameType);
        % Apply psychoacoustic model
        SMR = zeros(69,2);
        SMR(:,1) = psycho(frameT(:, 1), AACSeq3(k).frameType, frameTprevious1(:, 1), frameTprevious2(:, 1));
        SMR(:,2) = psycho(frameT(:, 2), AACSeq3(k).frameType, frameTprevious1(:, 2), frameTprevious2(:, 2));
        % Calculate T
        load('TableB219.mat');
        w_low = B219a(:,2);
        w_high = B219a(:,3);
        P = zeros(69,1);
        for b = 1:69 
            AACSeq3(k).chl.T(b) = sum(frameF((w_low(b)+1):(w_high(b)+1), 1).^2, 1);
            AACSeq3(k).chr.T(b) = sum(frameF((w_low(b)+1):(w_high(b)+1), 2).^2, 1);
        end
        AACSeq3(k).chl.T = AACSeq3(k).chl.T ./ SMR(:,1);
        AACSeq3(k).chr.T = AACSeq3(k).chr.T ./ SMR(:,2);
        % Apply quantizer
        S = zeros(1024,2);
        sfc = zeros(69,2);
        [S(:,1), sfc(:,1), AACSeq3(k).chl.G] = AACquantizer(frameF(:,1), AACSeq3(k).frameType, SMR(:,1));
        [S(:,2), sfc(:,2), AACSeq3(k).chr.G] = AACquantizer(frameF(:,2), AACSeq3(k).frameType, SMR(:,2));
        % Apply Huffman
        [AACSeq3(k).chl.stream, AACSeq3(k).chl.codebook] = encodeHuff(S(:,1), huffLUT);
        [AACSeq3(k).chr.stream, AACSeq3(k).chr.codebook] = encodeHuff(S(:,2), huffLUT);
        [AACSeq3(k).chl.sfc, ~] = encodeHuff(sfc(:,1), huffLUT);
        [AACSeq3(k).chr.sfc, ~] = encodeHuff(sfc(:,2), huffLUT); 
    end
    
end  

% Last printing and return
fprintf(repmat('\b',1,lineLength))
lineLength = fprintf('%2.0f / %4.0f',K-1,K-1);
fprintf('\n');
save(fnameAACoded, 'AACSeq3');


end