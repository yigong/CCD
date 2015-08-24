%EndOfTrack.m
load psft.mat
pixsize = 10.5;
dEdx = 0.22;
PSS = 0.25;
CSS = 0.25;
CW = 10;
%various tests of pixilation, interpolation, width and dE fitting


if true
    %pixilation along width of cuts
    
    x0 = (0:1:200)';
    y0 = zeros(size(x0));
    z0 = zeros(size(x0));
    dE0 = dEdx*ones(size(x0));
    
%     z = 0:50:650;
%     yoffset = 0:0.1:10.5;
    %quick:
    z = 0:650/4:650;
    yoffset = 0;
    D = cell(length(z),length(yoffset));
    cutX = cell(size(D));
    cutY = cell(size(D));
    cutE = cell(size(D));
    w = cell(size(D));
    dE = cell(size(D));
    
    progressbar(0);
    for i=1:length(z)
        for j=1:length(yoffset)
            %diffuse
            D{i,j} = DiffuseTrack4(psft,'manual',[x0,y0+yoffset(j),z0+z(i),dE0]);
            %make cuts
            cutx = mean(D{i,j}.cheat.xp(:,1)):PSS:(max(D{i,j}.cheat.xp(:,1))+CW/2);
            cuty = mean(D{i,j}.cheat.xp(:,2))-CW/2:CSS:mean(D{i,j}.cheat.xp(:,2))+CW/2;
            [cutX{i,j},cutY{i,j}] = meshgrid(cutx,cuty);
            cutE{i,j} = interp2(1:size(D{i,j}.img,1),1:size(D{i,j}.img,2),D{i,j}.img',cutX{i,j},cutY{i,j});
            %measure cuts
            for k=1:length(cutx)
                wfit = fit(cuty',cutE{i,j}(:,k),'gauss1',...
                    'StartPoint', [max(cutE{i,j}(:,k)),mean(cuty),CW/5], ...
                    'Lower', [0.5*max(cutE{i,j}(:,k))-.01,mean(cuty)-CW/2,.25], ...
                    'Upper', [2*max(cutE{i,j}(:,k))+.01,mean(cuty)+CW/2,CW/2]);
                w{i,j}(k) = 1.665 * wfit.c1 * pixsize;
                dE{i,j}(k) = sum(cutE{i,j}(:,k)) * CSS / pixsize; %not really dE, but dE/dx in keV/um.
            end
            progressbar(((i-1)*length(yoffset)+j)/(length(z)*length(yoffset)));
        end
    end
    
end