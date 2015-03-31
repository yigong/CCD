function [str]=rmquote(str)
%function removes "'s and replaces them w/ char39, which is '.  Useful for eval()
for j=1:length(str)
    if str(j)=='"'
        str(j)=char(39);
    end
end

