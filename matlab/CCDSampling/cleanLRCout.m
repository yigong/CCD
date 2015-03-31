for i_pos = 0:15
    NT = 0;
    alpha = [];
    alphaTrue = 0;
    beta = [];
    energy = [];
    ridgeWidth = {};
    wrongEnd = [];
    for iFile = 1:3
        fname = sprintf('./TrackerOuts%d_%d.mat', i_pos, iFile);
        if (i_pos == 9) & (iFile == 3)
            continue
        end
        
        try
            tmp = load(fname);
        catch err
            if err.identifier == 'MATLAB:load:couldNotReadFile'
                continue
            else
                'unknown error'
            end
        end
        NT = tmp.NT;
        alpha = [alpha;tmp.alpha];
        alphaTrue = tmp.alphaTrue;
        beta = [beta; tmp.beta];
        energy = [energy; tmp.energyMeasured];
        ridgeWidth = [ridgeWidth; tmp.ridgeWidth];
        wrongEnd = [wrongEnd; tmp.wrongEnd];
    end
    save(sprintf('./pos_%d', i_pos), 'NT', 'alpha', 'alphaTrue', 'beta', 'energy', 'ridgeWidth', 'wrongEnd')
end
