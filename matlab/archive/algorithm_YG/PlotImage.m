function PlotImage(img, cmap)
% plot function
font = 18;
figPos = [851, 388, 560, 420];
h = SurfElectronTrack(img, 'cmap', cmap);
set(h, 'fontsize', font);
set(get(h, 'parent'), 'position', figPos);
drawnow;




