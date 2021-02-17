function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)
% function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)
% AACQUANTIZER calculates hearing threshold T(b) and implements
% quantization for one channel
%
% Inputs:
% frameF      : Frame in frequency domain [128 x 8] for ESH, 
%               [1024 x 1] else
% frameType   : Type of frame according to SSC.m [String]
% SMR         : Signal to mask ratio from psychoacoustic model
%              [42 x 8] for ESH, [69 x 1] else
%
% Outputs:
% S           : [1024 x 1] with the quantized symbols of frame in Frequency
%               domain
% sfc         : Scale factor coefficients [42 x 8] for ESH, [69 x 1] else
% G           : Global gain of frame [8 x 1] for ESH, [1 x 1] else

load('TableB219.mat');
if frameType == "ESH"
    w_low = B219b(:,2); % [42 x 1]
    w_high = B219b(:,3); % [42 x 1]
    Nb = 42;
    S = zeros(128,8); % [128 x 8]->[1024 x 1]
    sfc = zeros(42,8); % [42 x 8]
    G = zeros(8,1); % [8 x 1]
        
    for i = 1:8
        X = frameF(:,i);
        % Calculate acoustic threshold 
        P = zeros(Nb, 1);
        for b = 1:Nb 
            P(b) = sum(X((w_low(b)+1):(w_high(b)+1)).^2);
        end
        T = P ./ SMR(:,i); % [69 x 1] 
        T(isnan(T)) = 0;
        % Find scale factor gain a
        a = zeros(Nb,1);
        for b = 1:Nb
            a(b) = round((16/3)*log2((max(X)^(3/4))/8191));
        end
        Pe = zeros(Nb, 1);
        idx = 1:Nb;
        while(true)
            if max(abs(a(2:end) - a(1:end-1))) > 60
                break; 
            end
            for b = idx % for each band calculate Pe(b) from S(k), X_hat(k)                
                S([(w_low(b)+1):(w_high(b)+1)],i) = sign(X([(w_low(b)+1):(w_high(b)+1)])).*fix((abs(X([(w_low(b)+1):(w_high(b)+1)])).*2.^(-1/4*a(b))).^(3/4)+0.4054);
                X_hat = sign(S([(w_low(b)+1):(w_high(b)+1)],i)).*(abs(S([(w_low(b)+1):(w_high(b)+1)],i)).^(4/3)).*2.^(1/4*a(b));              
                Pe(b) = sum((X([(w_low(b)+1):(w_high(b)+1)]) - X_hat).^2);
            end  
            idx = (find(Pe < T))';
            if length(idx) == 0
                break; 
            end
            a(idx) = a(idx) + 1;
        end
        sfc(:,i) = a;
        G(i) = a(1);
        for b = 2:Nb
            sfc(b,i) = a(b)-a(b-1);
        end
    end
    S = reshape(S, [1024 1]);
    S(isinf(S)|isnan(S)) = 0;
    G(isinf(G)|isnan(G)) = 0;
    sfc(isinf(sfc)|isnan(sfc)) = 0;
    return
else
    w_low = B219a(:,2); % [69 x 1]
    w_high = B219a(:,3); % [69 x 1]
    Nb = 69;
    S = zeros(1024,1); % [1024 x 1] 
    sfc = zeros(69,1); % [69 x 1]
    G = zeros(1,1); % [1 x 1]
    
    X = frameF;
    % Calculate acoustic threshold 
    P = zeros(Nb, 1);
    for b = 1:Nb 
        P(b) = sum(X((w_low(b)+1):(w_high(b)+1), :).^2);
    end
    T = P ./ SMR; % [69 x 1]
    T(isnan(T)) = 0;
    % Find scale factor gain a
    a = zeros(Nb,1);
    for b = 1:Nb
        a(b) = round((16/3)*log2((max(X)^(3/4))/8191));
    end
    Pe = zeros(Nb, 1);
    idx = 1:Nb;
    X_hat = zeros(size(X));
    while(true)
        if max(abs(a(2:end) - a(1:end-1))) > 60
            break; 
        end
        X_hat = zeros(size(X));
        for b = idx % for each band calculate Pe(b) from S(k), X_hat(k)
            Pe(b) = 0;
            for k = (w_low(b)+1):(w_high(b)+1) % re-calculate S(k), X_hat(k)
                S(k) = sign(X(k)).*fix((abs(X(k)).*2.^(-1/4*a(b))).^(3/4)+0.4054);
                X_hat(k) = sign(S(k)).*(abs(S(k)).^(4/3)).*2.^((1/4)*a(b));
            end
            Pe(b) = sum((X([(w_low(b)+1):(w_high(b)+1)]) - X_hat([(w_low(b)+1):(w_high(b)+1)])).^2); 
        end  
        idx = (find(Pe < T))'; % remove saturated bands
        if length(idx) == 0
            break; 
        end
        a(idx) = a(idx) + 1; % for the bands remaining, increase a(b) by 1
    end
    sfc = a;
    G = a(1);
    % DPCM
    for b = 2:Nb 
        sfc(b) = a(b)-a(b-1);
    end
    S(isinf(S)|isnan(S)) = 0;
    G(isinf(G)|isnan(G)) = 0;
    sfc(isinf(sfc)|isnan(sfc)) = 0;
    return;

end

