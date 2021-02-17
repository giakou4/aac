function [frameF] = filterbank(frameT, frameType, winType)
% function [frameF] = filterbank(frameT, frameType, winType)
% FILTERBANK implements MDCT after multiplication of the frame in time
% domain with a window.
%
% Inputs:
% frameT     : Frame in time domain [2048 x 2]
% frameType  : Type of frame according to SSC.m [String]
% winType    : Window type: KBD|SIN [String]
%
% Outputs:
% frameF     : Frame in frequency domain
%              [1024 x 2] for OLS|LSS|LPPS or
%              [128 x 2 x 8] for ESH
    
if frameType == "OLS"
   
    % Create Window
    if winType == "KBD"
        W = kbdwin(2048,6) 
    elseif winType == "SIN"
        W =  [sin(pi/2048*([0:1023]+0.5))';sin(pi/2048*([1024:2047]+0.5))'];
    else
        error('Error. Invalid Window Type. Choose among KBD|SIN');
    end
    
    % Take MDCT of the inner product of window and frame in time
    frameF = zeros(1024, 2);
    frameF = reshape(mdct(frameT,W,'PadInput',false),[1024 2]);  % from [1024x1x2] to [1024x2]
    return
    
elseif frameType == "LSS"
    
    % Create Window
    if winType == "KBD"
        WL = kbdwin(2048,6); WL = WL(1:1024);
        WR = kbdwin(256,4);  WR = WR(129:end);
        W = [WL; ones(448, 1); WR; zeros(448, 1)];
    elseif winType == "SIN"
        W = [sin(pi/2048*([0:1023]+0.5))'; ones(448, 1); sin(pi/256*([(128:255)]+0.5))'; zeros(448, 1)];
    else
        error('Error. Invalid Window Type. Choose among KBD|SIN');
    end
   
    % Take MDCT of the inner product of window and frame in time
    frameF = zeros(1024, 2);
    frameF = reshape(mdct(frameT,W,'PadInput',false),[1024 2]);  % from [1024x1x2] to [1024x2]
    return
    
elseif frameType == "ESH"
    
    % Create Window
    if winType == "KBD"
        W = kbdwin(256,4); % MATLAB 2020b and afterwards :(
    elseif winType == "SIN"
        W = [sin(pi/256*([0:127]+0.5))';sin(pi/256*([128:255]+0.5))'];
    else
        error('Error. Invalid Window Type. Choose among KBD|SIN');
    end
    
    % Initialize return: special case 8 matrices
    frameF = zeros(128,2,8); 
    for k=1:8
            % Take MDCT of the inner product of window and frame in time
            frameF(:,:,k) = reshape(mdct(frameT((448+128*k-127):(448+128*k+128),:),W,'PadInput',false),[128 2]);
    end
    return
    
elseif frameType == "LPS"
    
   % Create Window
    if winType == "KBD"
        WL = kbdwin(256,4);  WL = WL(1:128);
        WR = kbdwin(2048,6); WR = WR(1025:end);
        W = [zeros(448, 1); WL; ones(448, 1); WR];
    elseif winType == "SIN"
        W = [zeros(448, 1); sin(pi/256*([0:127]+0.5))'; ones(448, 1); sin(pi/2048*([(1024:2047)]+0.5))'];
    else
        error('Error. Invalid Window Type. Choose among KBD|SIN');
    end
   
    % Take MDCT of the inner product of window and frame in time
    frameF = zeros(1024, 2);
    frameF = reshape(mdct(frameT,W,'PadInput',false),[1024 2]);  % from [1024x1x2] to [1024x2]
    return 
    
end

end