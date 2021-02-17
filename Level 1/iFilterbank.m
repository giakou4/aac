function [frameT] = iFilterbank(frameF, frameType, winType)
% function [frameT] = iFilterbank(frameF, frameType, winType)
% IFILTERBANK implements the inverse procedure of Filterbank.
%
% Inputs:
% frameF     : Frame in frequency domain 
%              [1024 x 2] for OLS|LSS|LPPS or
%              [128 x 2 x 8] for ESH for ESH
% frameType  : Type of frame according to SSC.m [String]
% winType    : Window type: KBD|SIN [String]
%
% Outputs:
% frameT     : Frame in time domain [2048 x 2]

frameT = zeros(2048, 2);

if frameType == "OLS"
    
    % Create Window
    if winType == "KBD"
        W = kbdwin(2048,6) 
    elseif winType == "SIN"
        W =  [sin(pi/2048*([0:1023]+0.5))';sin(pi/2048*([1024:2047]+0.5))'];
    else
        error('Error. Invalid Window Type. Choose among KBD|SIN');
    end
    
    % Take IMDCT of frame in frequency and then multiply with window
    frameT(:, 1) = imdct(frameF(:,1),W,'PadInput',false); 
    frameT(:, 2) = imdct(frameF(:,2),W,'PadInput',false); 
    return
    
elseif frameType == "LSS"
    
    % Create Window
    if winType == "KBD"
        w = kbdwin(2048,6); WL = w(1:1024);
        w = kbdwin(256,4);  WR = w(129:end);
        W = [WL; ones(448, 1); WR; zeros(448, 1)];
    elseif winType == "SIN"
        W = [sin(pi/2048*([0:1023]+0.5))'; ones(448, 1); sin(pi/256*([(128:255)]+0.5))'; zeros(448, 1)];
    else
        error('Error. Invalid Window Type. Choose among KBD|SIN');
    end
    
    % Take IMDCT of frame in frequency and then multiply with window
    frameT(:, 1) = imdct(frameF(:,1),W,'PadInput',false); 
    frameT(:, 2) = imdct(frameF(:,2),W,'PadInput',false); 
    return
    
elseif frameType == "ESH"
    
    % Create Window
    if winType == "KBD"
        W = kbdwin(256,4); 
    elseif winType == "SIN"
        W = [sin(pi/256*([0:127]+0.5))';sin(pi/256*([128:255]+0.5))'];
    else
        error('Error. Invalid Window Type. Choose among KBD|SIN');
    end
    
    for channel=1:2
        for k=1:8
            frameT((((448+128*k-127):(448+128*k+128))), channel) = ...
                frameT((((448+128*k-127):(448+128*k+128))), channel) + ...
                imdct(frameF(:,channel,k),W,'PadInput',false);
        end
    end
    return
    
elseif frameType == "LPS"

    % Create Window
    if winType == "KBD"
        w = kbdwin(2048,6); WR = w(1025:end);
        w = kbdwin(256,4);  WL = w(1:128);
        W = [zeros(448, 1); WL; ones(448, 1); WR];
    elseif winType == "SIN"
        W = [zeros(448, 1); sin(pi/256*([0:127]+0.5))'; ones(448, 1); sin(pi/2048*([(1024:2047)]+0.5))'];
    else
        error('Error. Invalid Window Type. Choose among KBD|SIN');
    end
    
    % Take IMDCT of frame in frequency and then multiply with window
    frameT(:, 1) = imdct(frameF(:,1),W,'PadInput',false); 
    frameT(:, 2) = imdct(frameF(:,2),W,'PadInput',false); 
    return
    
end

end