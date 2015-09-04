function res = PeakFitting(b,n,ROI)
%-------------------------------------------------------------------------%
% ATTN SSRDM Students:
% For the b and n inputs use the output variables of the "hist" function.
% The syntax is:
%
% [n,x] = hist(energies_vector,[array_of_bin_centers]);
% bar(x,n)
%
% x = bin centers
% n = counts in each bin
%
% the array of bin centers should allow for your energy spectrum to have
% fine enough binning such that your photo-peaks have a shape and aren't a
% couple of bars on the bar plot. The bar plot is just a histogram.
%
% The output of the function: res = FHWM/peak_channel. This is the relative
% FWHM.
%-------------------------------------------------------------------------%

% simplifying additions for SSRDM students J.Curtis
% 07/28/2010
numPeaks = 1;
sig_est = 1;
plotflag = 0;
cubicflag = 0;

% This function was assembled by Dan Chivers.

% res = PeakFit(b,n,ROI,numPeaks,sig_est,plotflag)
%
%   General Peak Fitting Tool using Curve Fitting Toolbox.  
%       Single Peak Fitting Function: (Linear Background)
%       m*x + b + A/sqrt(2*pi*sig^2)*exp(-(x-pk)^2/(2*sig^2))
%
%       Dual Peak Fitting Function: (Linear Background)
%       m*x + b + A/sqrt(2*pi*sig^2)*exp(-(x-pk)^2/(2*sig^2))+A2/sqrt(2*pi*sig^2)*exp(-(x-pk2)^2/(2*sig^2))
%
%       Inputs:
%       b: Bin Centers (x-axis ... usually energy) vector
%       n: Counts vector
%       ROI: 2 element vector region of interest
%           data includes ROI values. (ie. b >= ROI(1) & b <= ROI(2))
%       numPeaks: 1 for single, 2 for double
%       sig_est: estimate for standard deviation of a single peak
%       plotflag: 0 (no Display) >0 (show display)
%
%   res: Results structure
%         res.ft = fit results
%         res.gof = goodness of fit results
%         res.ci = confidence limits matrix (1 sigma)
%         res.cnts1 = total counts in first peak
%         res.cnts2 = total counts in second peak (if applicable)
%         res.bkgcnts = total background counts
%         res.src1cnts = total source counts in peak 1
%         res.src2cnts = total source counts in peak 2 (if applicable)
%
%   3/19/2010 - Script created from SDV peak fitting base - DHC
%
%   3/22/2010 - Added Cubic background model



    if length(ROI) ~= 2
        error('ROI must have two elements [bmin,bmax]');
        return;
    end
    
    if size(b,1) > 1
        b = b';
    end
    
    if size(n,1) > 1
        n = n';
    end
    
    if nargin < 7
        cubicflag = 0;
    end

    b_o = b;
    n_o = n;
    
    F = find(b >= ROI(1) & b <= ROI(2));
    b = b(F); n = n(F);
    
    if numPeaks==1
        mxb = b(find(n==max(n)));
        mxb = mxb(1);

        %Estimates
        %sig_est = 35;  %keV
        
        A_est = max(n)*sqrt(2*pi*sig_est^2);
        
        if length(n) > 6
            m_est = (mean(n(1:3)) - mean(n(end-2:end)))/(b(2) - b(end-1));
        else
            m_est = 0;
        end
        b_est = b(1) - m_est*b(1); 

        s = fitoptions('Method','NonlinearLeastSquares',...
                            'Lower',[0,-inf,-inf,0,0],...
                            'Upper',[inf,inf,inf,inf,inf],...
                            'Startpoint',[A_est,b_est,m_est,mxb,sig_est]);
        f = fittype('m*x + b + A/sqrt(2*pi*sig^2)*exp(-(x-pk)^2/(2*sig^2))','options',s);
        [ft,gof] = fit(b',n',f);
        
        
        %Addind cubic fitting here
        % parameters: [A,b,m,m2,m3,pk,sig]
        
        if cubicflag
            s = fitoptions('Method','NonlinearLeastSquares',...
                                'Lower',[0,-inf,-inf,-inf,-inf,0,0],...
                                'Upper',[inf,inf,inf,inf,inf,inf,inf],...
                                'Startpoint',[ft.A,ft.b,ft.m,0,0,ft.pk,ft.sig]);
            f = fittype('m3*x^3 + m2*x^2 + m*x + b + A/sqrt(2*pi*sig^2)*exp(-(x-pk)^2/(2*sig^2))','options',s);
            [ft,gof] = fit(b',n',f);

        end   
        
        if plotflag
            b_hires = min(b):(b(2)-b(1))/10:max(b);
            semilogy(b_o,n_o,'b',b_o,n_o,'k.'); hold on;
            if ~cubicflag
                semilogy(b_hires,ft.b + ft.m*b_hires,'r','LineWidth',2.0);
                semilogy(b_hires,ft.b + ft.m*b_hires + ...
                    ft.A/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk).^2/(2*ft.sig^2)),'g','LineWidth',2.0);
            else
                semilogy(b_hires,ft.b + ft.m*b_hires + ft.m2*b_hires.^2 + ft.m3*b_hires.^3,'r','LineWidth',2.0);
                semilogy(b_hires,ft.b + ft.m*b_hires + ...
                    ft.m2*b_hires.^2 + ft.m3*b_hires.^3 + ...
                    ft.A/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk).^2/(2*ft.sig^2)),'g','LineWidth',2.0);
            end
            drawnow;
            hold off;
        end
        
        ci = confint(ft,.68);
        cnts = sum(n);
        if cubicflag
            bkg_ = ft.b + ft.m*b + ft.m2*b.^2 + ft.m3*b.^3;
            os = 2;
        else
            bkg_ = ft.b + ft.m*b;
            os = 0;
        end
        bkg = round(sum(bkg_));
        src = round(sum(n - bkg_));

        disp('Fit Results');
        disp('-----------');
        disp(['Peak Value: ',num2str(ft.pk),' ( ',num2str(ci(1,4+os)),' , ',num2str(ci(2,4+os)),' )']);
        disp(['FWHM: ',num2str(2.355*ft.sig),' ( ',num2str(2.355*ci(1,5+os)),' , ',num2str(2.355*ci(2,5+os)),' )']);
        disp(['Total Counts: ',num2str(cnts),' ( ',num2str(round(cnts-sqrt(cnts))),' , ',num2str(round(cnts+sqrt(cnts))),' )']);
        disp(['Background Counts: ',num2str(bkg),' ( ',num2str(round(bkg-sqrt(bkg))),' , ',num2str(round(bkg+sqrt(bkg))),' )']);
        disp(['Peak Counts: ',num2str(src),' ( ',num2str(round(src-sqrt(src+bkg))),' , ',num2str(round(src+sqrt(src+bkg))),' )']);
        disp(['Peak Fit Area: ',num2str(ft.A/(b(2)-b(1))),' ( ',num2str(ci(1,1)/(b(2)-b(1))),' , ',num2str(ci(2,1)/(b(2)-b(1))),' )']);

        res = 2.355*ft.sig/ft.pk;
        disp(['Relative FWHM = ',num2str(2.355*ft.sig./ft.pk.*100),'%']);
        
        res.ft = ft;
        res.ci = ci;
        res.gof = gof;
        res.cnts1 = cnts;
        res.bkgcnts = bkg;
        res.src1cnts = src;
        
        
    elseif numPeaks==2
        
        m_est = (n(end)-n(1))/(b(end)-b(1));
        b_est = n(1) - m_est*b(1);

        NBkg = b_est + m_est*b;
        DBkg = n - NBkg;
        
        F = find(DBkg < 0);
        if ~isempty(F)
            indMax = find(DBkg == min(DBkg));
            indMax = indMax(1);
            R = -1*DBkg(indMax)/n(indMax);
            m_est = (1-R)*m_est;
            b_est = (1-R)*b_est;
        end

        ind1 = find(DBkg == max(DBkg));
        A1_est = n(ind1) - (b_est + m_est*b(ind1));
        sig1_est = sig_est;
        n1 = b_est + m_est*b + A1_est*exp(-(b-b(ind1)).^2/(2*sig1_est^2));
        n2 = n-n1;
        ind2 = find(n2 == max(n2));
        A2_est = n(ind2) - (b_est + m_est*b(ind2));
      



        s = fitoptions('Method','NonlinearLeastSquares',...
                        'Lower',[0,0,-inf,-inf,0,0,0],...
                        'Upper',[inf,inf,inf,inf,inf,inf,inf],...
                        'Startpoint',[A1_est,...
                                      A2_est,...
                                      b_est,...
                                      m_est,...
                                      b(ind1),...
                                      b(ind2),...
                                      sig1_est,...
                                      ]);
        f = fittype('m*x + b + A/sqrt(2*pi*sig^2)*exp(-(x-pk)^2/(2*sig^2))+A2/sqrt(2*pi*sig^2)*exp(-(x-pk2)^2/(2*sig^2))','options',s);
        [ft,gof] = fit(b',n',f);
        
        
        if cubicflag

            s = fitoptions('Method','NonlinearLeastSquares',...
                        'Lower',[0,0,-inf,-inf,-inf,-inf,0,0,0],...
                        'Upper',[inf,inf,inf,inf,inf,inf,inf,inf,inf],...
                        'Startpoint',[ft.A,...
                                      ft.A2,...
                                      ft.b,...
                                      ft.m,...
                                      0,...
                                      0,...
                                      ft.pk,...
                                      ft.pk2,...
                                      ft.sig,...
                                      ]);
            f = fittype('m3*x^3 + m2*x^2 + m*x + b + A/sqrt(2*pi*sig^2)*exp(-(x-pk)^2/(2*sig^2))+A2/sqrt(2*pi*sig^2)*exp(-(x-pk2)^2/(2*sig^2))','options',s);
            [ft,gof] = fit(b',n',f);      
        end
        
        
        if plotflag
            b_hires = min(b):(b(2)-b(1))/10:max(b);
            semilogy(b_o,n_o,'b',b_o,n_o,'k.'); hold on
            if ~cubicflag
                semilogy(b_hires,ft.b + ft.m*b_hires,'r','LineWidth',2.0);
                semilogy(b_hires,ft.b + ft.m*b_hires + ft.A/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk).^2/(2*ft.sig^2)),'m','LineWidth',2.0);
                semilogy(b_hires,ft.b + ft.m*b_hires + ft.A2/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk2).^2/(2*ft.sig^2)),'g','LineWidth',2.0);
                semilogy(b_hires,ft.b + ft.m*b_hires + ft.A/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk).^2/(2*ft.sig^2))+ft.A2/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk2).^2/(2*ft.sig^2)),'k','LineWidth',2.0);
            else
                semilogy(b_hires,ft.b + ft.m*b_hires + ft.m2*b_hires.^2 + ft.m3*b_hires.^3,'r','LineWidth',2.0);
                semilogy(b_hires,ft.b + ft.m*b_hires + ft.m2*b_hires.^2 + ft.m3*b_hires.^3 + ft.A/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk).^2/(2*ft.sig^2)),'m','LineWidth',2.0);
                semilogy(b_hires,ft.b + ft.m*b_hires + ft.m2*b_hires.^2 + ft.m*b_hires.^3 + ft.A2/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk2).^2/(2*ft.sig^2)),'g','LineWidth',2.0);
                semilogy(b_hires,ft.b + ft.m*b_hires + ft.m2*b_hires.^2 + ft.m3*b_hires.^3 + ft.A/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk).^2/(2*ft.sig^2))+ft.A2/sqrt(2*pi*ft.sig^2)*exp(-(b_hires-ft.pk2).^2/(2*ft.sig^2)),'k','LineWidth',2.0);
            end    
            hold off
            drawnow;
        end
        
        cnts = sum(n);
        ci = confint(ft,.68);
        if ~cubicflag
            n1 = ft.b + ft.m*b + ft.A/sqrt(2*pi*ft.sig)*exp(-(b-ft.pk).^2/(2*ft.sig^2));
            n2 = ft.b + ft.m*b + ft.A2/sqrt(2*pi*ft.sig)*exp(-(b-ft.pk2).^2/(2*ft.sig^2));
        else
           n1 = ft.b + ft.m*b + ft.m2*b.^2 + ft.m3*b.^3 + ft.A/sqrt(2*pi*ft.sig)*exp(-(b-ft.pk).^2/(2*ft.sig^2));
            n2 = ft.b + ft.m*b + ft.A2/sqrt(2*pi*ft.sig)*exp(-(b-ft.pk2).^2/(2*ft.sig^2));
        end
            
        db = b(2)-b(1);
        
        cnts1 = sum(n1)/db;
        cnts2 = sum(n2)/db;
        if ~cubicflag
            bkg_ = ft.b + ft.m*b;
            os = 0;
        else
            bkg_ = ft.b + ft.m*b + ft.m2*b.^2 + ft.m3*b.^3;
            os = 2;
        end
        
        bkg = round(sum(bkg_));
        src1 = round(sum(n1 - bkg_));
        src2 = round(sum(n2 - bkg_));

        disp('Fit Results');
        disp(['Background Counts: ',num2str(bkg),' ( ',num2str(round(bkg-sqrt(bkg))),' , ',num2str(round(bkg+sqrt(bkg))),' )']);
        disp(['Peak 1 Value: ',num2str(ft.pk),' ( ',num2str(ci(1,5+os)),' , ',num2str(ci(2,5+os)),' )']);
        disp(['FWHM 1: ',num2str(2.355*ft.sig),' ( ',num2str(2.355*ci(1,7+os)),' , ',num2str(2.355*ci(2,7+os)),' )']);
        disp(['Peak 1 Fit Area: ',num2str(ft.A/db),' ( ',num2str(ci(1,1)/db),' , ',num2str(ci(2,1)/db),' )']);
        disp(['Peak 2 Value: ',num2str(ft.pk2),' ( ',num2str(ci(1,6)),' , ',num2str(ci(2,6)),' )']);
        disp(['FWHM 2: ',num2str(2.355*ft.sig),' ( ',num2str(2.355*ci(1,7+os)),' , ',num2str(2.355*ci(2,7+os)),' )']);
        disp(['Peak 2 Fit Area: ',num2str(ft.A2/db),' ( ',num2str(ci(1,2)/db),' , ',num2str(ci(2,2)/db),' )']);
        disp(['Total Counts: ',num2str(cnts),' (',num2str(cnts-sqrt(cnts)),' , ',num2str(cnts+sqrt(cnts)),' )']);

        res = 2.355*ft.sig/ft.pk;
        disp(['Relative FWHM = ',num2str(2.355*ft.sig./ft.pk.*100),'%']);
        
%         JCurtis 07/30/2010         
%         res.ft = ft;
%         res.gof = gof;
%         res.ci = ci;
%         res.cnts1 = cnts1;
%         res.cnts2 = cnts2;
%         res.bkgcnts = bkg;
%         res.src1cnts = src1;
%         res.src2cnts = src2;
%         res.cubicflag = cubicflag;
        
    end
    