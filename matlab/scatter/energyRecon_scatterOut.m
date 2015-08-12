figure();
set(gcf, 'Position', [440, 378, 640, 400])

Er = cell(1, 3);
lc = containers.Map(thk, {'b','r', 'g'});
mc2 = 511.0;
for i = 1:3
    Er{i} = mc2./(cosPhi{i}.*(sqrt(1+2*mc2./E{i})) - 1);
    [counts, bins] = hist(Er{i}, [1649:1:1801]);
    plot(bins(2:end-1), counts(2:end-1)/max(counts(2:end-1)), ...
        'Color', lc(thk{i}), 'LineWidth', 2)
    hold on 
end
set(gca,'FontSize',18)
xlabel('E_r(keV) 1 keV width', 'FontSize', 18)
ylabel('Counts', 'FontSize', 18)
legend('10um', '30um', '50um')
grid on