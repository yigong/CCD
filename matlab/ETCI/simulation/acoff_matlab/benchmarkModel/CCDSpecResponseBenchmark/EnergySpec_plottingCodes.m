    % Calculate Initial Number of Gamma's
    nps = 101*500000; %4*10*10^6+2*5*10^6; %Total Number of Starting Photons
 
    
    figure(2)
    subplot(2,1,2);
    plot(y1, x1./nps, 'k', ...
         y1, x2./nps, 'b', ...
         y1, x3./nps, 'r','LineWidth',2);

    legend('Total Energy Spectrum','Cointained Electrons','Non-contained Electrons');
    title('Simulated CCD Spectrum of Experiment', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]','Fontsize',14);
    ylabel('Event Probability','Fontsize',14);
    %axes('Fontsize',20)
    nps2 = length(E);
    subplot(2,1,1);plot(y1,x4./nps2,'b','LineWidth',2);
    title('Experimental CCD Total Energy Spectrum','Fontsize',14);  
    xlabel('CCD Event Energy [keV]','Fontsize',14);
    ylabel('Event Probability [arb units]','Fontsize',14);