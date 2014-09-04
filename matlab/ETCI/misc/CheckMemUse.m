function b = CheckMemUse(w)
%function b = CheckMemUse(w)
%
% get total memory usage in bytes
%
% Input: w = whos;
%

%sum up bytes from w
b = 0;
for i=1:length(w)
    b = b + w(i).bytes;
end
