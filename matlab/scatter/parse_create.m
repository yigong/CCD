matrl = cell(1, 3);
matrl{1} = 'Li_';
matrl{2} = 'Nitro_';
matrl{3} = 'Pb_';
cosPhi = cell(1, 3);
E = cell(1, 3);

for j = 2
    for i = 0:14
        fname = ['/Users/Yigong/Google Drive/Research/CCD/matlab/scatter/',...
            matrl{j}, num2str(i), '.mat'];
        tmp = load(fname);
        cosPhi{j} = [cosPhi{j}, tmp.cosPhi_create];
        E{j}      = [E{j}, tmp.E_create];
    end
end