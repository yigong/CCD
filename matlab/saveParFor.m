function saveParFor(fName, var, ver)
% save to use in parfor to avoid transparency violation
    save(fName, 'var', ver);
end
