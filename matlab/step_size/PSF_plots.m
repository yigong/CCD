dataDir = '/Users/Yigong/Google Drive/Research/CCD/matlab/step_size/data';
cd(dataDir);
fileTmp = dir('PSF_10nm_backP.mat')'
figure()
set(gcf, 'Position', [440, 378, 640, 480])

j = 0
keys = {'backP', 'pixelP'};
values = {'-', '--'};
ls = containers.Map(keys, values); % linestyle
keys = {'50um', '10um', '3um', '1um', '300nm', '100nm', '60nm', '30nm', '10nm'};
values = {[1, 0.27, 0], [1, 0.4, 0.6], 'r', 'm', 'g', 'c', 'b', [0.41, 0.41, 0.41], 'b'};
lc = containers.Map(keys, values); % line color

legendIdx= [];    % legend index
legendStr = cell(size(fileTmp));    % legend strings
h = [];    % plot handles
for i = 1:length(fileTmp)
    fileName = fileTmp(i).name;
    splitted = strsplit(fileName, {'_','.'});
    step = splitted{2};
    plane = splitted{3};
    load(fileName);
    % f_pos = sqrt((xT - xM).^2 + (yT - yM).^2) < 50;
    %f_E   = abs(ET - EM) < 1;
    f_Emin = ET > 0;
    % f_al  = (alphaM < 90) | (alphaM > 270);
    f = f_Emin;
    [counts, bins] = hist(xWindow(f)/1000, [-4.2:0.2:4.2]);
    counts = counts(2:end-1);
    counts = counts;
    bins = bins(2:end-1);
    
    h(i) = plot(bins+j, counts, 'Color', lc(step), 'LineStyle', ls(plane), 'LineWidth', 2);
    legendStr{i} = step;
    if strcmp(plane, 'backP')
        legendIdx = [legendIdx, i];
        
        %set(h{i}, 'DisplayName', step)
    end
    
    hold on
    clearvars alphaM alphaT betaM EM ET xM xT xWindow yM yT f
end
%legend(h(legendIdx), legendStr(legendIdx))
xlabel('Photon Position (mm)', 'FontSize', 18)
ylabel('Counts', 'FontSize', 18)
ylim([0, 730])
%line([-0.2, -0.2], [0, 730], 'Color', 'r', 'LineWidth', 2)
%line([0.2, 0.2], [0, 730], 'Color', 'r', 'LineWidth', 2)

%l1 = legend(gca, 'show')
set(gca, 'FontSize',18)
grid on