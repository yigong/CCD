thk = cell(1, 3);
thk{1} = 'thk10um_';
thk{2} = 'thk30um_';
thk{3} = 'thk50um_';
cosPhi = cell(1, 3);
E = cell(1, 3);
x = cell(1, 3);
y = cell(1, 3);
z = cell(1, 3);

for j = 1:3
    for i = 0:14
        fname = ['/Users/Yigong/Google Drive/Research/CCD/matlab/scatter/',...
            thk{j}, num2str(i), '.mat'];
        tmp = load(fname);
        cosPhi{j} = [cosPhi{j}, tmp.cosPhi_scatterOut ];
        E{j}      = [E{j}, tmp.E_scatterOut];
        x{j}      = [x{j}, tmp.x_scatterOut];
        y{j}      = [y{j}, tmp.y_scatterOut];
        z{j}      = [z{j}, tmp.z_scatterOut - 3500];
    end
end
