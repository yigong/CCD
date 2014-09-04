function M2 = processmuon(M2,MC,tdt)
%function M2 = processmuon(M2,MC,tdt)
%
% M2: cell array of muon structures
% MC: MuonCoord
% tdt: total drift time table (TotalDriftTime.mat)
%
%Muon structure:
% .DT = drift time. Lookup from tdt.
% .var = variance = (w/2.355)^2
% .DC = diffusion constant fit from variance & drift time
% .DCrmse = RMSE on diffusion constant fit
% .yint = y-intercept of DC fit
%
%NOT IMPLEMENTED:
% .z2 = Extrapolate diffusion constant down to one pixel (19.9 um^2), recalculate depths.
% .DT2 = from z2, lookup in tdt.
% .DC2 = diffusion constant again.

tstep = 0.01;   %time step for tdt table
progressbar(0); 
for i=1:length(M2); 
    if ~isempty(M2{i}); 
        M2{i}.DT = interp1(0:tstep:650,tdt(:,2),M2{i}.z); 
        %add variance = (FWHM / 2.355)^2
        M2{i}.var = (M2{i}.w *10.5/2.355).^2;
        
        %fit diffusion constant: var = 2*D*(DT)
        
        x = length(M2{i}.var);
        %19.9 = variance for a FWHM of one pixel = minimum variance we
            %expect to see
        o = fitoptions('poly1','upper',[Inf,19.9],'lower',[0,19.9]);
        %skip MC(i,1)==0
        if ~MC(i,1)==0
            xdata = M2{i}.DT(logical(~isnan(M2{i}.w(1:x))));
            ydata = M2{i}.var(logical(~isnan(M2{i}.w(1:x))));
            [f,gof] = fit(xdata',ydata','poly1',o);
            M2{i}.DC = f.p1 / 2;
%             M2{i}.yint = f.p2;
            M2{i}.DCrmse = gof.rmse;
%             M2{i}.z2 = 
%         elseif MC(i,1)==-1
%             f = fit(M2{i}.DT(x:-1:1)',M2{i}.var(x:-1:1)','poly1',o);
%             M2{i}.DC = f.p1 / 2;
        end
        
    end
    
    progressbar(i/length(M2)); 
end