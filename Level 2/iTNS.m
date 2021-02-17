function frameFout = iTNS(frameFin, frameType, TNScoeffs)
% function [frameFout, TNScoeffs] = TNS(frameFin, frameType)
% TNS implements Temporal Noise Shaping for 1 channel
%
% Inputs
%
% frameFin     : MDCT coefficients [128 x 8] for "ESH", [1024 x 1] else
%                after TNS 
% frameType    : Type of frameFin
% TNScoeffs    : Quantized TNS coefficients [4 x 8] for "ESH", [4 x 1] else
%
% Outputs
% frameFout    : MDCT coefficients [128 x 8] for "ESH", [1024 x 1] else
%                before TNS

frameFout = zeros(size(frameFin));
if frameType == "ESH"
    for i = 1:8 % for each subframe
        isStable = isstable([1, -TNScoeffs(:, i)']);
        if isStable == 0
            warning('Filtre is not stable');
        end
        frameFout(:, i) = filter(1, [1, -TNScoeffs(:, i)'], frameFin(:, i));
    end
    return
else
    isStable = isstable([1, -TNScoeffs']);
    if isStable == 0
        warning('Filtre is not stable');
    end
    frameFout = filter(1, [1, -TNScoeffs'], frameFin);
end
    

end

