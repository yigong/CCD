function dAlpha = FixAlphaError(dAlpha)
%function dAlpha = FixAlphaError(dAlpha)
%
% Shifts values onto the range of [-180, 180].
%
% dAlpha input can be any numerical matrix.
% dAlpha output is the same size as input.

if ~isnumeric(dAlpha)
    error('dAlpha input should be numeric.')
end

valuesHaveChanged = true;
while valuesHaveChanged
    oldDAlpha = dAlpha;
    % Fix finite values only
    outOfRange = abs(dAlpha)>180 & isfinite(dAlpha);
    adjustSign = dAlpha(outOfRange) ./ abs(dAlpha(outOfRange));
    dAlpha(outOfRange) = dAlpha(outOfRange) - 360*adjustSign;
    
    valuesHaveChanged = ...
        ~all(dAlpha(isfinite(dAlpha))==oldDAlpha(isfinite(oldDAlpha)));
end