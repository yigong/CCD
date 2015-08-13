figure();
set(gcf, 'Position', [440, 378, 640, 400])
Er = cell(1, 3);
lc = containers.Map(matrl, {'r', 'b', 'g'});
mc2 = 511.0;
for i = 2
    Er{i} = mc2./(cosPhi{i}.*(sqrt(1+2*mc2./E{i})) - 1);
    [counts, bins] = hist(Er{i}(1:39000), [1728.9:0.1:1735.1]);
    plot(bins(2:end-1), counts(2:end-1), ...
        'Color', lc(matrl{i}), 'LineWidth', 2)
    hold on 
end
set(gca,'FontSize',18)
xlabel('E_r(keV)', 'FontSize', 18)
ylabel('Counts', 'FontSize', 18)
%legend('Nitrocellulose Z=7')
grid on