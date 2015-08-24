function saveParFor(fName, result, ver)
% save to use in parfor to avoid transparency violation
    save(fName, '-struct', 'result', ver);
end
