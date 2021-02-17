function SNR = demoAAC1(fNameIn,fNameOut)
% function SNR = demoAAC1(fNameIn, fNameOut)
% DEMOAAC2 demonstrates the procedure in Level 2
%
% Inputs
% fnameIN    : Name of file in .wav format. Sound must have 2 channels with
%              sampling frequency equals to 48kHz
% fNameOut   : Name of the output .wav file. Sound must have 2 channels 
%              with sampling frequency equals to 48kHz
%
% Outputs
% SNR        : Signal to Noise Ration in dB of procedure in Level 2 
%
% Example:
% demoAAC1('LicorDeCalandraca.wav', 'LicorDeCalandraca_demoAAC1.wav')
tic

[x,Fs] = audioread(fNameIn);

AACSeq1 = AACoder1(fNameIn);
xdec = iAACoder1(AACSeq1, fNameOut);

% Do the same as AACoder1 
N = length(x);
extraSamples = 1024 - mod(N, 1024);
x = [zeros(2048,2); x; zeros(2048+extraSamples,2)];
x = x(1:end-1024,:);

SNR = zeros(2,1);
for channel = 1:2
    SNR(channel) = snr(x(:,channel),xdec(:,channel)-x(:,channel));
end

toc
%{
if ~exist([pwd, '\..\'],'dir' )
    mkdir('Report')
end
dir = [pwd, '\..\Report\'];
for channel = 1:2
    figure
    plot(x(:,channel));
    hold on;
    plot(xdec(:,channel));
    title(['Channel ',num2str(channel)]); 
    xlabel('samples');
    ylabel('Amplitude'); 
    legend('initial signal','decoded signal');
    saveas(gcf,[dir 'demoAAC1-ch',num2str(channel),'.png']) 

    figure
    plot(x(:,channel)-xdec(:,channel));
    title(['Error of channel ',num2str(channel)]); 
    saveas(gcf,[dir 'demoAAC1-error-ch',num2str(channel),'.png']) 
end
%}
end

