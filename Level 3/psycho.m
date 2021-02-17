function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)
% function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)
% PSYCHO implemets psychoacoustic model for 1 channel.
%
% Inputs:
% frameT      : Frame in time domain [2048 x 2]
% frameType   : Type of frame according to SSC.m [String]
% frameTprev1 : Previous frame [2048 x 2]
% frameTprev2 : Previous of previous frame [2048 x 2]
%
% Outputs
% SMR         : Signal to mask ration [42 x 8] for ESH, [69 x 1] else

% Load whatever is necessary
load('TableB219.mat');
if frameType == 'ESH'
    w_low = B219b(:,2); % [42 x 1]
    w_high = B219b(:,3); % [42 x 1]
    bval = B219b(:,5); % [42 x 1]
    qsthr = B219b(:, 6); % [42 x 1]
    Nb = 42;
    spreadingfunction = zeros(Nb,Nb);
    SMR = zeros(42,8);
else
    w_low = B219a(:,2); % [69 x 1]
    w_high = B219a(:,3); % [69 x 1]
    bval = B219a(:,5); % [69 x 1]
    qsthr = B219a(:, 6); % [69 x 1]
    Nb = 69;
    spreadingfunction = zeros(Nb,Nb);
    SMR = zeros(69,1);
end

% Spreading function, calculate again and again ... :(
for i = 1:Nb
    for j = 1:Nb
        if i >= j
            tmpx = 3*(bval(j)-bval(i));
        else
            tmpx = 1.5*(bval(j)-bval(i));
        end
        tmpz = 8*min(0,(tmpx-0.5)^2-2*(tmpx-0.5));
        tmpy = 15.811389 + 7.5*(tmpx+0.474) - 17.5*(1+(tmpx+0.474)^2)^(0.5);
        if tmpy < -100
            spreadingfunction(i,j) = 0;
        else
            spreadingfunction(i,j) = 10^((tmpz+tmpy)/10);
        end
    end
end

% Check type 
if frameType == 'ESH'
    frames = zeros(256,10); % [256 x 10]
    frames(:, 1) = frameTprev1(1217:1472); % We just need the previous subframe - note that we ignore right 448
    frames(:, 2) = frameTprev1(1345:1600); % and previous of previous subframe - note that we ignore right 448: 448+1600=2048
    for k = 1:8
        frames(:, k+2) = frameT(448+(((k-1)*128+1):((k+1)*128)));
    end
    for k = 3:10 % for each of the 8 subframes, k is the current subframe
        % Calculate c from r, r_pred and f_pred
        r = zeros(128,3); % r(:,1) for current, r(:,2:3) for prev1,2
        f = zeros(128,3); % f(:,1) for current, f(:,2:3) for prev1,2
        for i = 0:2
            s = frames(:,k-i); % i=0 -> current
            sw = s.*(0.5 - 0.5*cos(pi/256*([1:256]-0.5)))';
            Sw = fft(sw);
            r(:,i+1) = abs(Sw(1:128));
            f(:,i+1) = angle(Sw(1:128));
        end
        r_pred = 2*r(:,2) - r(:,3);
        f_pred = 2*f(:,2) - f(:,3);
        r = r(:,1);
        f = f(:,1);
        C = sqrt((r.*cos(f) - r_pred.*cos(f_pred)).^2 ...
             + (r.*sin(f) - r_pred.*sin(f_pred)).^2)./(r + abs(r_pred));    
        
        % Calculate  energy and predictability for each b
        e = zeros(Nb,1);
        c = zeros(Nb,1);
        for b = 1:Nb
            e(b) = sum(r((w_low(b)+1):(w_high(b)+1),:).^2);
            c(b) = sum(C((w_low(b)+1):(w_high(b)+1),:).*r((w_low(b)+1):(w_high(b)+1),:).^2);
        end

        % Combine energy and predictability with spreading function function
        ecb = (e'*spreadingfunction)';
        ct = (c'*spreadingfunction)';

        % Normalize
        cb = ct ./ ecb;
        en = zeros(Nb,1);
        for b = 1:Nb
            en(b) = ecb(b) / (sum(spreadingfunction(:,b)));
        end

        % Calculate tonality index in [0,1]
        tb = -0.299 - 0.43 * log(cb);
        tb(tb>1) = 1;
        tb(tb<0) = 0;
        %tb = (tb - min(tb))/(max(tb)-min(tb));

        % Calculate SNR
        SNR = tb*18 + 6*(1-tb);

        % Transform from dB to energy ratio
        bc = 10.^(-SNR/10);

        % Calculate energy threshold
        nb = en .* bc;

        % Calculate noise level
        qthr_est = eps()*(256/2)*10.^(qsthr/10);
        npart = max([nb,qthr_est]')';

        % Calculate Signal to Mask Ratio
        SMR(:,k-2) = e ./ npart;         
    end
    return

else
    % calculate c from r, r_pred and f_pred
    frames = [frameT,frameTprev1,frameTprev2]; % [2048 x 3] current, prev, prev-prev
    r = zeros(1024,3); % [1024 x 3]
    f = zeros(1024,3); % [1024 x 3]
    for i=1:3 % for each frameT of selected channel
        s =  frames(:,i);
        sw = s.*(0.5 - 0.5*cos(pi/2048*([1:2048]-0.5)))';
        Sw = fft(sw);
        r(:,i) = abs(Sw(1:1024));
        f(:,i) = angle(Sw(1:1024));
    end
    r_pred = 2*r(:,2) - r(:,3);
    f_pred = 2*f(:,2) - f(:,3);
    r = r(:,1);
    f = f(:,1);
    C = sqrt((r.*cos(f)-r_pred.*cos(f_pred)).^2 ...
         + (r.*sin(f)-r_pred.*sin(f_pred)).^2)./(r+abs(r_pred));
    
    % Calculate  energy and predictability for each b
    e = zeros(Nb,1);
    c = zeros(Nb,1);
    for b = 1:Nb
        e(b) = sum(r((w_low(b)+1):(w_high(b)+1),:).^2);
        c(b) = sum(C((w_low(b)+1):(w_high(b)+1),:).*r((w_low(b)+1):(w_high(b)+1),:).^2);
    end
    
    % Combine energy and predictability with spreading function
    ecb = (e'*spreadingfunction)';
    ct = (c'*spreadingfunction)';
    
    % Normalize
    cb = ct ./ ecb;
    en = zeros(Nb,1);
    for b = 1:Nb
        en(b) = ecb(b) / (sum(spreadingfunction(:,b)));
    end
    
    % Calculate tonality index in [0,1]
    tb = -0.299 - 0.43 * log(cb);
    tb(tb>1) = 1;
    tb(tb<0) = 0;
    %tb = (tb - min(tb))/(max(tb)-min(tb));
    
    % Calculate SNR
    SNR = tb*18 + 6*(1-tb);
    
    % Transform from dB to energy ratio
    bc = 10.^(-SNR/10);
    
    % Calculate energy threshold
    nb = en .* bc;
    
    % Calculate noise level
    qthr_est = eps()*(2048/2)*10.^(qsthr/10);
    npart = max([nb,qthr_est]')';
    
    % Calculate Signal to Mask Ratio
    SMR = e ./ npart;
    return;
end

end

