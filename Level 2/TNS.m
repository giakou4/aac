function [frameFout, TNScoeffs] = TNS(frameFin, frameType)
% function [frameFout, TNScoeffs] = TNS(frameFin, frameType)
% TNS implements Temporal Noise Shaping for 1 channel
%
% Inputs
%
% frameFin     : MDCT coefficients [128 x 8] for "ESH", [1024 x 1] else
% frameType    : Type of frameFin
%
% Outputs
% frameFout    : MDCT coefficients [128 x 8] for "ESH", [1024 x 1] else
%                after TNS
% TNScoeffs    : Quantized TNS coefficients [4 x 8] for "ESH", [4 x 1] else

load('TableB219.mat');
if frameType == "ESH"
    b = B219b(:,2); %w_low from w2203tfa.pdf TableΒ.2.1.9.b 48kHz short FFT
    b(end+1)=128;
    Nb = length(b);
    frameFout = zeros(128,8);
    TNScoeffs = zeros(4,8);
    for i = 1:8 % for each subframe
        % Use X as temp and moreover keep the symbols from assignment :)
        X = frameFin(:, i);
        % Energy for each band j
        P = zeros(Nb-1,1);
        for j = 1:Nb-1
            P(j) = sum(X((b(j)+1):b(j+1)).^2);
        end
        % Normalization factor Sw
        Sw = zeros(128, 1);
        for j = 1:Nb-1
            Sw((b(j)+1):b(j+1)) = sqrt(P(j));
        end
        % Smoothing
        for k = (length(Sw)-1):-1:1
            Sw(k) = (Sw(k) + Sw(k+1))/2; 
        end
        for k = 2:(length(Sw))
            Sw(k) = (Sw(k) + Sw(k-1))/2; 
        end
        % Normalize MDCT coefficients
        Xw = X./Sw;
        % Compute a from solving R*a=r 
        [~, H] = corrmtx(Xw, 4);
        r = H(2 : 5, 1);
        R = H(1 : 4, 1 : 4);
        a = (R \ r)';
        % Quantize a
        step = 0.1;
        bits = 4;
        range = step*(2^bits - 1);
        a_min = -range/2;
        a_max = range/2;
        a(a > a_max) = a_max;
        a(a < a_min) = a_min; 
        a = round(a/(step))*step;
        % Get TNS coefficients
        TNScoeffs(:, i) = a';
        % Check if filter is stable
        isStable = isstable([1, -TNScoeffs(:, i)']);
        if isStable == 0
            warning('Filtre is not stable');
        end
        % Get new frame in frequency 
        frameFout(:, i) = filter([1, -TNScoeffs(:, i)'], 1, X);
    end
else
    b = B219a(:,2); %w_low from w2203tfa.pdf TableΒ.2.1.9.a 48kHz long FFT
    b(end+1) = 1024;
    Nb = length(b);
    frameFout = zeros(1024,1);
    TNScoeffs = zeros(4,1);
    X = frameFin;
    % Energy for each band j
    P = zeros(Nb-1,1);
    for j = 1:(Nb-1)
        P(j) = sum(X((b(j)+1):b(j+1)).^2);
    end
    % Normalization factor Sw
    Sw = zeros(1024, 1);
    for j = 1:(Nb-1)
        Sw((b(j)+1):b(j+1)) = sqrt(P(j));
    end
    % Smoothing
    for k = (length(Sw)-1):-1:1
        Sw(k) = (Sw(k) + Sw(k+1))/2; 
    end
    for k = 2:1:(length(Sw))
        Sw(k) = (Sw(k) + Sw(k-1))/2; 
    end
    % Normalize MDCT coefficients
    Xw = X./Sw;  
    Xw(isnan(Xw)) = 0;
    % Compute a from solving R*a=r 
    [~, H] = corrmtx(Xw, 4);
    r = H(2 : 5, 1);
    R = H(1 : 4, 1 : 4);
    a = (R\r)';
    a(isnan(a)) = 0;
    % Quantize a
    step = 0.1;
    bits = 4;
    range = step*(2^bits - 1);
    a_min = -range/2;
    a_max = range/2;
    a(a > a_max) = a_max;
    a(a < a_min) = a_min; 
    a = round(a/(step))*step;
    % Get TNS coefficients
    TNScoeffs = a';
    % Check if filter is stable
    isStable = isstable([1, -TNScoeffs']);
    if isStable == 0
        warning('Filtre is not stable');
    end
    % Get new frame in frequency 
    frameFout = filter([1, -TNScoeffs'], 1, X);
end

end

