function SNR = demoAAC2(fNameIn,fNameOut)
% function SNR = demoAAC2(fNameIn, fNameOut)
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
% demoAAC2('LicorDeCalandraca.wav', 'LicorDeCalandraca_demoAAC2.wav')
tic

[x,Fs] = audioread(fNameIn);

AACSeq2 = AACoder2(fNameIn);
xdec = iAACoder2(AACSeq2, fNameOut);

% Do the same as AACoder2 
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
    saveas(gcf,[dir 'demoAAC2-ch',num2str(channel),'.png']) 

    figure
    plot(x(:,channel)-xdec(:,channel));
    title(['Error of channel ',num2str(channel)]); 
    saveas(gcf,[dir 'demoAAC2-error-ch',num2str(channel),'.png']) 
end
%}
end

