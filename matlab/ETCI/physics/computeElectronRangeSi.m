function rangeUm = computeElectronRangeSi(energyKev,reference)
% function rangeUm = computeElectronRangeSi(energyKev,reference)
% 
% Compute the extrapolated range in silicon, for electrons of the given energy.
% 
% Input:
%   energyKev: electron energy in KeV. This can be of any dimension.
%   reference: (optional) 1996, 2002, or 'CSDA'. 
%       1996 and 2002 refer to extrapolated range, calculated from the formulae 
%           in Tabata's two (more recent) papers, and using Tabata's formula for
%           CSDA range too. See below for discussion. 
%       'CSDA' returns the CSDA range, from Tabata's formula from 1996a.
%       (Default: 2002) (1996 and 2002 can be numeric or string inputs)
% 
% Output:
%   rangeUm: range in um, in an array of size(energyKev).
% 
% Discussion of Tabata references
%   First of all, the CSDA range calculation is good from 10 keV to 100 MeV. It 
% is validated to EPSTAR and ESTAR from NIST and IAEA, and the error is 
% negligible compared to the extrapolated range expressions. From 1 to 10 keV, 
% the CSDA range is uncertain due to shell effects which are not accounted for.
% Error relative to EPSTAR is 0.7% RMS and 1.5% maximum. (see [1])
%   Tabata 1996a [1] fits a formula using Monte Carlo data from ITS from 
% 100 keV to 100 MeV for elements from Be to U, and extended with experimental
% data from 10 keV to 100 keV. Error relative to IPS is 0.9% RMS, 2.0% max (max 
% occurs at Al at 100 keV). Error relative to experiment is 14% RMS, 28% max. 
% ([1], Table 3)
%   Tabata 2002 [2] fits an identical formula, for elements Be, C, Al, Cu, Ag, 
% Au, U with PENELOPE data between 100 keV and 50 MeV, and extends it with data 
% published by various authors from 1923 to 1964 (see sec 2.3). The RMS error 
% from PENELOPE is 1.5% (over the range from 100 keV to 50 MeV), and maximum 
% error from PENELOPE for Al is 1.5% at 100 keV. The RMS error from experiment
% for 10 keV to 100 keV is 18% because of discrepancies among different
% measurements. It seems that Al is one of the better fits though.
%   1996a range is below 2002 range by less than 3.3% above 250 keV, but by as 
% much as as 12% at 15 keV. (Above 25 MeV 1996a range is higher.)
%   In conclusion, 1996a might be better above 50 MeV (where I don't care) or
% below 10 keV (where CSDA range is bad anyway). 2002 is fitted with equivalent
% methods, and more recently with better data.
% 
%   Bottom line: Use 2002.
% 
% [1] T. Tabata et al., Nucl. Instr. and Meth. in Phys. Res. B 119 
%       (1996) 463--470
% [2] T. Tabata et al., Radiation Physics and Chemistry 64 (2002) 161--167

% let computeElectronRange.m handle input

% silicon constants
Z = 14;             
A = 28.085;         %amu
I_EV = 171;         %eV
RHO_GCM3 = 2.329;    % g/cm^3

rangeUm = computeElectronRange(energyKev,Z,A,I_EV,RHO_GCM3,reference);