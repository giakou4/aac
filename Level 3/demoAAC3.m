function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, frameAACoded)
% function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, frameAACoded)
% DEMOAAC3 demonstrate the procedure in Level 3
%
% Inputs:
% fNameIn      : Name of file in .wav format. Sound must have 2 channels with
%                sampling frequency equals to 48kHz
% fnameOut     : Name of the output .wav file. Sound must have 2 channels 
%                with sampling frequency equals to 48kHz
% fnameAACoded : Name of .mat file to save AACSeq3
%
% Outputs:
% SNR          : Signal to Noise Ration in dB of procedure in Level 3 
% bitrate      : bits per second
% compression  : bitrate before compression divided with bitrate after.
%
% Example:
% [SNR, bitrate, compression] = demoAAC3('LicorDeCalandraca.wav','LicorDeCalandraca_demoAAC3.wav','LicorDeCalandraca_AACSeq3_demoAAC3.mat')

tic
[x,Fs] = audioread(fNameIn);

AACSeq3 = AACoder3(fNameIn,frameAACoded);
xdec = iAACoder3(AACSeq3, fNameOut);

% Do the same as AACoder2 
N = length(x);
extraSamples = 1024 - mod(N, 1024);
x = [zeros(2048,2); x; zeros(2048+extraSamples,2)];
x = x(1:end-1024,:);

% Calculate SNR
SNR = zeros(2,1);
for channel = 1:2
    SNR(channel) = snr(x(:,channel),xdec(:,channel)-x(:,channel));
end

% Load .mat file to remove T in order to compute its size
load(frameAACoded);
for i = 1:length(AACSeq3)
    AACSeq3(i).chl.T = [];
    AACSeq3(i).chr.T = [];
end

% Save New .mat file without T
save(['New_' frameAACoded], 'AACSeq3');

% Calculate bitrate before and after compression and the compresion ratio
frameAACoded_dir = dir(['New_' frameAACoded]);
frameAACoded_Size = frameAACoded_dir.bytes*8;
bitrate_after_compression = (frameAACoded_Size*Fs)/N;
bitrate = bitrate_after_compression;

fNameIn_dir = dir(fNameIn);
fNameIn_Size = fNameIn_dir.bytes*8;
bitrate_before_compression = (fNameIn_Size*Fs)/N;

compression = fNameIn_Size/frameAACoded_Size;
toc

%{
if ~exist([pwd, '\..\'],'dir' )
    mkdir('Report')
end
dir = [pwd, '\..\Report\'];
for channel = 1:2
    figure
    plot(xdec(:,channel));
    hold on;
    plot(x(:,channel));
    title(['Channel ',num2str(channel)]); 
    xlabel('samples');
    ylabel('Amplitude'); 
    legend('initial signal','decoded signal');
    saveas(gcf,[dir 'demoAAC3-ch',num2str(channel),'.png']) 

    figure
    plot(x(:,channel)-xdec(:,channel));
    title(['Error of channel ',num2str(channel)]); 
    saveas(gcf,[dir 'demoAAC3-error-ch',num2str(channel),'.png']) 
end
%}
end

