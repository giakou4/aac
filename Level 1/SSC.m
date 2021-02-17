function frameType = SSC(frameT, nextFrameT, prevFrameType)
% function frameType = SSC(frameT, nextFrameT, prevFrameType)
% FRAMETYPE implements Sequence Segmentation Control
%
% Inputs:
% frameT        : Frame in time domain. Includes 2 channels [2048 x 2]
% nextFrameT    : Next frame in time domain. Includes 2 channels [2048 x 2]
% prevFrameType : Type of previous frame [String]
%
% Output:
% frameType     : Type of frame [String]
%                 "OLS": for ONLY_LONG_SEQUENCE
%                 "LSS": for LONG_START_SEQUENCE
%                 "ESH": for EIGHT_SHORT_SEQUENCE
%                 "LPS": for LONG_STOP_SEQUENCE

if prevFrameType == "LPS"
    frameType = "OLS";
    return
elseif prevFrameType == "LSS"
    frameType = "ESH";
    return
end

for channel = 1:2
    
    % Calculate Energy and Attack Values of Next Filtered Frame
    filteredNextFrameT = filter([0.7548, -0.7548], [1, -0.5095], nextFrameT(:, channel));
    energy = zeros(8,1);
    for l=1:8
        energy(l) = sum(filteredNextFrameT(((l-1)*256+1):(l*256)).^2);
    end
    attackvalues = zeros(8, 1); 
    for l = 2:8
        attackvalues(l) = l*energy(l)/sum(energy(1:(l-1))); % attackvalue(1) = 0
    end
        
    % Find if Next Frame is ESH or not
    if sum(energy>10^(-3))>0 && sum(attackvalues>10)
        isESH = 1;
    else
        isESH = 0;
    end
    
    % Save option of frame type
    if prevFrameType == "OLS"
        if(isESH == 1)
            possibleFrameType(channel) = "LSS";
        else
            possibleFrameType(channel) = "OLS";
        end
    elseif prevFrameType == "ESH"
        if(isESH == 1)
            possibleFrameType(channel) = "ESH";
        else
            possibleFrameType(channel) = "LPS";
        end
    else % for initialization where prevFrameType = ""
        if (isESH == 1)
            possibleFrameType(channel) = "LSS";
        else
            possibleFrameType(channel) = "OLS";
        end
    end
    
end
   
% Sum up results for each channel and make a decision
if possibleFrameType(1) == "OLS" & possibleFrameType(2) == "OLS"
    frameType = "OLS";
    return
elseif (possibleFrameType(1) == "OLS" & possibleFrameType(2) == "LSS") || (possibleFrameType(1) == "LSS" & possibleFrameType(2) == "OLS") || (possibleFrameType(1) == "LSS" & possibleFrameType(2) == "LSS") 
    frameType = "LSS";
    return
elseif (possibleFrameType(1) == "OLS" & possibleFrameType(2) == "LPS") || (possibleFrameType(1) == "LPS" & possibleFrameType(2) == "OLS") || (possibleFrameType(1) == "LPS" & possibleFrameType(2) == "LPS")
    frameType = "LPS";
    return
else % all other...
    frameType = "ESH";
    return
end

end


