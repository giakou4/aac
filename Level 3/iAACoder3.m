function x = iAACoder3(AACSeq3, fNameOut)
% function x = iAACoder3(AACSeq3, fNameOut)
% IAACODER2 does the inverse job if AACoder2
%
% Inputs:
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
%                                       ['101010...']
%                * AACSeq3(i).chr.sfc - Huffman coded sequence of sfc
%                                       ['101010...']
%                * AACSeq3(i).chl.stream - Huffman coded sequence of frame  
%                                          in frequency domain ['101010']
%                * AACSeq3(i).chr.stream - Huffman coded sequence of frame
%                                          in frequency domain ['101010']
%                * AACSeq3(i).chl.codebook - Huffman's codebook [int]
%                * AACSeq3(i).chr.codebook - Huffman's codebook [int]
% fNameOut   : Name of the output .wav file. Sound must have 2 channels 
%              with sampling frequency equals to 48kHz
%
% Outputs:
% x          : if nargaout==1 contains the decoded sequence of samples
%
% Example:
% x = iAACoder3(AACSeq3, 'LicorDeCalandraca_demoAAC3.wav')

% Add path of Level 1 and 2
addpath([pwd '\..\Level 1'])
addpath([pwd '\..\Level 2'])

% Match AACoder3
warning off;
K = length(AACSeq3);
Fs = 48000;
winType = "SIN";
N = (K+1)*1024;
x = zeros(N,2);
huffLUT = loadLUT();
 
% Some printing
fprintf('Iteration: '); 
lineLength = fprintf('%2.0f / %4.0f',1,K);

% And so it begins
for k = 1:K
    if mod(k,10)==0
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf('%2.0f / %4.0f',k,K);
    end
    if AACSeq3(k).frameType == "ESH"
        % Get sfc by Huffman decoding
        sfc = zeros(42,8,2);
        sfcL = decodeHuff(AACSeq3(k).chl.sfc, 11, huffLUT); % 11 by observation
        sfc(:,:,1) = reshape(sfcL(1:42*8),[42 8]);
        sfcR = decodeHuff(AACSeq3(k).chr.sfc, 11, huffLUT); % 11 by observation
        sfc(:,:,2) =  reshape(sfcR(1:42*8),[42 8]);
        % Get S by Huffman decoding
        S = zeros(1024,2);
        SL = decodeHuff(AACSeq3(k).chl.stream, AACSeq3(k).chl.codebook, huffLUT);
        S(:,1) = reshape(SL(1:1024),[1024 1]);
        SR = decodeHuff(AACSeq3(k).chr.stream, AACSeq3(k).chr.codebook, huffLUT);
        S(:,2) = reshape(SR(1:1024),[1024 1]);
         % Get frameF after TNS by applying iAACquantizer
        frameF = zeros(128,8,2);
        frameF(:,:,1) = iAACquantizer(S(:,1), sfc(:,:,1), AACSeq3(k).chl.G, AACSeq3(k).frameType);
        frameF(:,:,2) = iAACquantizer(S(:,2), sfc(:,:,2), AACSeq3(k).chr.G, AACSeq3(k).frameType);
        % Get frameF before TNS by applying iTNS
        frameF(:,:,1) = iTNS(frameF(:,:,1),AACSeq3(k).frameType, AACSeq3(k).chl.TNScoeffs); 
        frameF(:,:,2) = iTNS(frameF(:,:,2),AACSeq3(k).frameType, AACSeq3(k).chr.TNScoeffs); 
        frameF = reshape(frameF,[128 2 8]); % could have done it [128 x 8 x 2] 
        % Get frameT by applying iFilterbank
        x(((k-1)*1024 + 1):(k+1)*1024, :) = x(((k-1)*1024 + 1):(k+1)*1024, :)+iFilterbank(frameF, AACSeq3(k).frameType, AACSeq3(k).winType);  
    else
        % Get sfc by Huffman decoding
        sfc = zeros(69,2);
        if k > 1
            sfcL = decodeHuff(AACSeq3(k).chl.sfc, 11, huffLUT); % 11 by observation
            sfc(:,1) = sfcL(1:69);
            sfcR = decodeHuff(AACSeq3(k).chr.sfc, 11, huffLUT); % 11 by observation
            sfc(:,2) = sfcR(1:69);
        end
        % Get S by Huffman decoding
        S = zeros(1024,2);
        if k > 1
            SL = decodeHuff(AACSeq3(k).chl.stream, AACSeq3(k).chl.codebook, huffLUT);
            S(:,1) = reshape(SL(1:1024),[1024 1]);
            SR = decodeHuff(AACSeq3(k).chr.stream, AACSeq3(k).chr.codebook, huffLUT);
            S(:,2) = reshape(SR(1:1024),[1024 1]);
        end
        % Get frameF after TNS by applying iAACquantizer
        frameF = zeros(1024,2);
        frameF(:,1) = iAACquantizer(S(:,1), sfc(:,1), AACSeq3(k).chl.G, AACSeq3(k).frameType);
        frameF(:,2) = iAACquantizer(S(:,2), sfc(:,2), AACSeq3(k).chr.G, AACSeq3(k).frameType);
        % Get frameF before TNS by applying iTNS
        frameF(:,1) = iTNS(frameF(:,1),AACSeq3(k).frameType, AACSeq3(k).chl.TNScoeffs); 
        frameF(:,2) = iTNS(frameF(:,2),AACSeq3(k).frameType, AACSeq3(k).chr.TNScoeffs); 
        % Get frameT by applying iFilterbank
        x(((k-1)*1024 + 1):(k+1)*1024,:) = x(((k-1)*1024 + 1):(k+1)*1024,:)+iFilterbank(frameF, AACSeq3(k).frameType, AACSeq3(k).winType);
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
audiowrite(fNameOut, x, Fs)

end

