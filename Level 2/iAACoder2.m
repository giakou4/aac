function x = iAACoder2(AACSeq2, fNameOut)
% function AACSeq2 = iAACoder2(fNameIn)
% IAACODER2 does the inverse job if AACoder2
%
% Inputs:
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
% fNameOut   : Name of the output .wav file. Sound must have 2 channels 
%              with sampling frequency equals to 48kHz
%
% Outputs:
% x          : if nargaout==1 contains the decoded sequence of samples
%
% Example:
% x = iAACoder2(AACSeq2, 'LicorDeCalandraca_demoAAC2.wav')

% Match AACoder2
warning off;
addpath([pwd '\..\Level 1']);
K = length(AACSeq2);
Fs = 48000;
winType = "SIN";
N = (K+1)*1024;
x = zeros(N,2);

% Some printing
fprintf('Iteration: '); 
lineLength = fprintf('%2.0f / %4.0f',1,K);

% And so it begins
for k = 1:K
    if mod(k,10)==0
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf('%2.0f / %4.0f',k,K);
    end
    if AACSeq2(k).frameType == "ESH"
        % Get frameF before TNS by applying iTNS
        frameF = zeros(128,2,8); % [128 x 2 x 8]
        frameF(:,1,:) = iTNS(AACSeq2(k).chl.frameF,...
                          AACSeq2(k).frameType, AACSeq2(k).chl.TNScoeffs);
        frameF(:,2,:) = iTNS(AACSeq2(k).chr.frameF,...
                          AACSeq2(k).frameType, AACSeq2(k).chr.TNScoeffs); 
        % Get frameT by applying iFilterbank
        frameT = iFilterbank(frameF, AACSeq2(k).frameType, winType);
            x((1024*k-1023):(1024*k+1024),:) = x((1024*k-1023):(1024*k+1024),:) + frameT; %overlapping is sum of right and left 
    else
        % Get frameF before TNS by applying iTNS
        frameF = zeros(1024,2); % [1024 x 2]
        frameF(:,1) = iTNS(AACSeq2(k).chl.frameF,AACSeq2(k).frameType, AACSeq2(k).chl.TNScoeffs);
        frameF(:,2) = iTNS(AACSeq2(k).chr.frameF,AACSeq2(k).frameType, AACSeq2(k).chr.TNScoeffs);
        % Get frameT by applying iFilterbank
        frameT = iFilterbank(frameF, AACSeq2(k).frameType, winType);
        x((1024*k-1023):(1024*k+1024),:) = x((1024*k-1023):(1024*k+1024),:) + frameT; %overlapping is sum of right and left 
  
    end
end

% Last printing
fprintf(repmat('\b',1,lineLength))
lineLength = fprintf('%2.0f / %4.0f',K,K);
fprintf('\n');

% Remove errors and write signal
x(x>1) = 1;
x(x<-1) = -1;
x(isnan(x)) = 0;
audiowrite(fNameOut, x, Fs);

end

