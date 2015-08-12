
% f_pos = sqrt((xT - xM).^2 + (yT - yM).^2) < 50;
% f_E   = abs(ET - EM) < 1;
f_Emin = ET > 500;
% f_al  = (alphaM < 90) | (alphaM > 270);
f = f_Emin & f_E;
f = logical(ones(size(ET)))
ss = 20; % scatter size
sty = 'filled'; % scatter style



figure();
subplot(2,3,1)
histogram(xWindow(f)/1000, [-5:0.1:5])
xlabel('mm')
title('PSF','FontSize',18)
set(gca,'FontSize',18)

subplot(2,3,2)
scatter(alphaT(f), alphaM(f), ss, sty)
xlabel('True (deg)','FontSize',18)
xlim([0, 360])
ylabel('Measured (deg)','FontSize',18)
ylim([0, 360])
title('Alpha','FontSize',18)
set(gca,'FontSize',18)

subplot(2,3,3)
scatter(ET(f), EM(f), ss , sty);
xlabel('True (keV)','FontSize',18)
xlim([0, 1700])
ylabel('Measured (keV)','FontSize',18)
title('Energy','FontSize',18)
set(gca,'FontSize',18)

subplot(2,3,4)
scatter(xT(f)/10.5, xM(f)/10.5, ss, sty);
xlabel('True (px)','FontSize',18)
ylabel('Measured (px)','FontSize',18)
title('X','FontSize',18)
set(gca,'FontSize',18)

subplot(2,3,5)
scatter(yT(f)/10.5, yM(f)/10.5, ss, sty);
xlabel('True (px)','FontSize',18)
ylabel('Measured (px)','FontSize',18)
title('Y','FontSize',18)
set(gca,'FontSize',18)

subplot(2,3,6)
scatter(xT(f)/10.5, yT(f)/10.5, ss, sty)
xlabel('True X (px)','FontSize',18)
ylabel('True Y (px)','FontSize',18)
title([num2str(sum(f)), ' events'],'FontSize',18)
set(gca,'FontSize',18)
