function x = iAACoder1(AACSeq1, fNameOut)
% function x = iAACoder1(AACSeq1, fNameOut) 
% IAACODER1 does the inverse procedure if AACoder1
%
% Inputs:
% AACSeq1    : Struct [K x 1] where K is the number of frames coded.
%              Contains:
%              * AACSeq1(i).frameType - OLS|LSS|LPPS|ESH [String]
%              * AACSeq1(i).winType - SIN|KBD [String]
%              * AACSeq1(i).chl.frameF - MDCT coefficient of left channel
%                                        [1024 x 1] for OLS|LSS|LPPS
%                                        [128 x 8] for ESH
%              * AACSeq1(i).chr.frameF - MDCT coefficient of right channel
%                                        [1024 x 1] for OLS|LSS|LPPS
%                                        [128 x 8] for ESH
% fNameOut    : Name of the output .wav file. Sound must have 2 channels 
%               with sampling frequency equals to 48kHz
%
% Outputs:
% x           : if nargaout==1 contains the decoded sequence of samples
% 
% Example:
% x = iAACoder1(AACSeq1, 'LicorDeCalandraca_demoAAC1.wav')

% Match AACoder1
warning off;
K = length(AACSeq1);
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
    if AACSeq1(k).frameType == "ESH"
        frameF = zeros(128,2,8); % [128 x 2 x 8]
        frameF(:,1,:) = AACSeq1(k).chl.frameF;
        frameF(:,2,:) = AACSeq1(k).chr.frameF;
        % Apply Inverse Filterbank to get Frame in Time
        frameT = iFilterbank(frameF, AACSeq1(k).frameType, winType);
        x((1024*k-1023):(1024*k+1024),:) = x((1024*k-1023):(1024*k+1024),:) + frameT; %overlapping is sum of right and left
    else
        frameF = zeros(1024,2); % [1024 x 2] 
        frameF(:,1) = AACSeq1(k).chl.frameF;
        frameF(:,2) = AACSeq1(k).chr.frameF;
        % Apply Inverse Filterbank to get Frame in Time
        frameT = iFilterbank(frameF, AACSeq1(k).frameType, AACSeq1(k).winType);
        x((1024*k-1023):(1024*k+1024),:) = x((1024*k-1023):(1024*k+1024),:) + frameT; %overlapping is sum of right and left
    end
end

% Last printing 
fprintf(repmat('\b',1,lineLength))
lineLength = fprintf('%2.0f / %4.0f',K,K);
fprintf('\n');

% Write signal
x(x>1) = 1;
x(x<-1) = -1;
x(isnan(x)) = 0;
audiowrite(fNameOut, x, Fs);
end
