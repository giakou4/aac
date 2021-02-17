function [frameF] = iAACquantizer(S, sfc, G, frameType)
% function [frameF] = iAACquantizer(S, sfc, G, frameType)
% iAACquantizer does the inverse job of AACquantizer
%
% Inputs:
% S           : [1024 x 1] with the quantized symbols of frame in Frequency
%               domain
% sfc         : Scale factor coefficients [42 x 8] for ESH, [69 x 1] else
% G           : Global gain of frame [8 x 1] for ESH, [1 x 1] else
% frameType   : Type of frame according to SSC.m
%
% Outsputs:
% frameF      : Frame in frequency domain [128 x 8] for "ESH", 
%               [1024 x 1] else

load('TableB219.mat');
if frameType == "ESH"
    w_low = B219b(:,2); % [42 x 1]
    w_high = B219b(:,3); % [42 x 1]
    Nb = 42;
    frameF = zeros(128,8); % [256 x 8]
    S = reshape(S, [128 8]);
    for i = 1:8
        a = zeros(Nb,1);
        % Inverse DPCM
        for b = 2:Nb
            a(b) = sfc(b,i) + a(b-1);
        end
        a = a + G(i);
        for b = 1:Nb
            for j = (w_low(b)+1):(w_high(b)+1)
                frameF(j, i) = sign(S(j, i)).*(abs(S(j,i)).^(4/3)).*2.^(1/4 * a(b));              
            end
        end
    end
else
    w_low = B219a(:,2); % [69 x 1]
    w_high = B219a(:,3); % [69 x 1]
    Nb = 69;
    frameF = zeros(1024,1); % [1024 x 1]
    a = zeros(Nb,1);
    % Inverse DPCM
    for b = 2:Nb
        a(b) = sfc(b) + a(b-1);
    end
    a = a + G;
    for b = 1:Nb 
        for j = (w_low(b)+1):(w_high(b)+1)
            frameF(j) = sign(S(j)).*(abs(S(j)).^(4/3)).*2.^(1/4 * a(b));    
        end
    end
end

end

