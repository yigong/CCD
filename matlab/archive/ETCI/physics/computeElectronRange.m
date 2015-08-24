function rangeUm = computeElectronRange(energyKev,Z,A,Iev,rhoGcm3,reference)
% function rangeUm = computeElectronRange(energyKev,Z,A,Iev,rhoGcm3,reference)
% 
% Compute the extrapolated range, for electrons of the given energy.
% 
% Input:
%   energyKev: electron energy in KeV. This can be of any dimension.
%   Z: atomic number or effective atomic number
%   A: average atomic weight, in amu
%   Iev: mean ionization energy. Si: 171, Ge: 350, CdTe: 539.3
%       **get this from http://physics.nist.gov/cgi-bin/Star/compos.pl
%   rhoGcm3: density in g/cm^3
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

% handle input
defaultReferenceMode = 2002;
if nargin==5
    reference = defaultReferenceMode;
elseif nargin<5
    error('Require five arguments...')
end
if (ischar(reference) && ...
        (strcmp(reference,'1996') || strcmp(reference,'1996a'))) || ...
        (isnumeric(reference) && reference==1996)
    referenceMode = '1996a';
elseif (ischar(reference) && strcmp(reference,'2002')) || ...
        (isnumeric(reference) && reference==2002)
    referenceMode = '2002';
elseif ischar(reference) && strcmpi(reference,'CSDA')
    referenceMode = 'csda';
else
    warning('Unrecognized reference argument... defaulting to 2002.')
    referenceMode = defaultReferenceMode;
end

% define what we use to convert CSDA to extrapolated range
deviationFactorFunction = assignFunction(referenceMode, energyKev);


% CSDA range. (Continuous Slowing-Down Approximation)
energyMev = energyKev ./ 1e3;
csdaRangeGCm2 = computeCsdaRange(energyMev, Z, A, Iev);
csdaRangeCm = csdaRangeGCm2 ./ rhoGcm3;
csdaRangeUm = csdaRangeCm .* 1e4;

% extrapolated range.
deviationFactor = deviationFactorFunction(energyMev, Z);
rangeUm = deviationFactor .* csdaRangeUm;     % Tabata 1996a equation 1


function functionHandle = assignFunction(referenceMode, energyKev)
% function functionHandle = assignFunction(referenceMode, energyKev)
switch referenceMode
    case '1996a'
        % formula validity constants
        MAX_ENERGY = 100e3;     % 100 MeV from Tabata 1996a p 2
        MIN_ENERGY = 1;         % 1 keV from Tabata 1996a p 2
        SEMI_MIN_ENERGY = 10;   % 10 keV, below which shell effects are relevant
                                % but not considered, from Tabata 1996a p 2
        
        % check for validity
        if any(energyKev < MIN_ENERGY)
            warning('Range is inaccurate below 1 keV.')
        end
        if any(energyKev > MIN_ENERGY & energyKev < SEMI_MIN_ENERGY)
            warning(['Range is less accurate below 10 keV, ' ...
                'because shell effects are ignored.'])
        end
        if any(energyKev > MAX_ENERGY)
            warning('Range is inaccurate above 100 MeV.')
        end
        
        functionHandle = @computeDeviationFactor1996a;
    case '2002'
        % formula validity constants
        MAX_ENERGY = 50e3;      % 50 MeV from Tabata 2002 sec 1; 2.1
        MIN_ENERGY = 10;        % 10 keV from Tabata 2002 sec 2.3
        
        % check for validity
        if any(energyKev < MIN_ENERGY)
            warning('Range is inaccurate below 10 keV. Try Tabata 1996a')
        end
        if any(energyKev > MAX_ENERGY)
            warning('Range is inaccurate above 50 MeV. Try Tabata 1996a')
        end
        
        functionHandle = @computeDeviationFactor2002;
    case 'csda'
        % formula validity constants
        MAX_ENERGY = 100e3;     % unclear but implied by 1996a
        MIN_ENERGY = 1;         % 1 keV from Tabata 1996a p 2
        SEMI_MIN_ENERGY = 10;   % 10 keV, below which shell effects are relevant
                                % but not considered, from Tabata 1996a p 2
        
        % check for validity
        if any(energyKev < MIN_ENERGY)
            warning('Range is inaccurate below 1 keV.')
        end
        if any(energyKev > MIN_ENERGY & energyKev < SEMI_MIN_ENERGY)
            warning(['Range is less accurate below 10 keV, ' ...
                'because shell effects are ignored.'])
        end
        if any(energyKev > MAX_ENERGY)
            warning('Range is inaccurate above 100 MeV.')
        end
        
        functionHandle = @unity;
end


function deviationFactor = computeDeviationFactor1996a(energyMev, Z)
% function deviationFactor = computeDeviationFactor1996a(energyMeV, Z)
% 

ELECTRON_REST_ENERGY_MEV = 0.510999;

t0 = energyMev ./ ELECTRON_REST_ENERGY_MEV;

% constants from Tabata's fit. Tabata 1996a Table 1
B = [...
    0.3879, ...     1
    0.2178, ...     2
    0.4541, ...     3
    0.03068, ...    4
    3.326e-16, ...  5
    13.24, ...      6
    1.316, ...      7
    14.03, ...      8
    0.7406, ...     9
    4.294e-3, ...   10
    1.684, ...      11
    0.2264, ...     12
    0.6127, ...     13
    0.1207, ...     14
    ];

% secondary values. Tabata 1996a equations 3--8
a = [...
    B(1) .* Z.^B(2), ...                    1
    B(3) + B(4).*Z, ...                     2
    B(5) .* Z.^(B(6) - B(7).*log(Z)), ...   3
    B(8) ./ Z.^(B(9)), ...                  4
    B(10).*Z.^(B(11) - B(12).*log(Z)), ...  5
    B(13) .* Z.^B(14), ...                  6
    ];

% Tabata 1996a equation 2:
deviationFactor = 1 ./ (a(1) + a(2)./(1 + a(3)./t0.^a(4) + a(5).*t0.^a(6)));


function deviationFactor = computeDeviationFactor2002(energyMev, Z)
% function deviationFactor = computeDeviationFactor2002(energyMeV, Z)
% 

ELECTRON_REST_ENERGY_MEV = 0.510999;

t0 = energyMev ./ ELECTRON_REST_ENERGY_MEV;

% constants from Tabata's fit. Tabata 2002 Table 2
B = [...
    0.2946, ...     1
    0.2740, ...     2
    18.4, ...       3
    3.457, ...      4
    1.377, ...      5
    6.59, ...       6
    2.414, ...      7
    1.094, ...      8
    0.05242, ...    9
    0.3709, ...     10
    0.958, ...      11
    2.02, ...       12
    1.099, ...      13
    0.2808, ...     14
    0.2042, ...     15
    ];

% secondary values. Tabata 2002 equations 9--14
a = [...
    B(1) .* Z.^B(2), ...                    1
    B(3) .* Z.^(-B(4) + B(5).*log(Z)), ...  2
    B(6) .* Z.^(-B(7) + B(8).*log(Z)), ...  3
    B(9) .* Z.^(B(10)), ...                 4
    B(11) .* Z.^(-B(12)+B(13).*log(Z)), ... 5
    B(14) .* Z.^B(15), ...                  6
    ];

% Tabata 2002 equation 8: (identical to 1996a eq. 2)
deviationFactor = 1 ./ (a(1) + a(2)./(1 + a(3)./t0.^a(4) + a(5).*t0.^a(6)));


function deviationFactor = unity(energyMev, ~)
% function deviationFactor = unity(energyMev, ~)
% 
% Return 1's, because CSDA range = 1 * CSDA range.

deviationFactor = ones(size(energyMev));


function csdaRangeGCm2 = computeCsdaRange(energyMev, Z, A, I_EV)
% function csdaRangeGCm2 = computeCsdaRange(energyMev, Z, A, IeV)
% 

ELECTRON_REST_ENERGY_MEV = 0.510999;

I = I_EV / (1e6*ELECTRON_REST_ENERGY_MEV);
t0 = energyMev ./ ELECTRON_REST_ENERGY_MEV;

% Tabata 1996a Table 2
D = [...
    3.600, ...      1
    0.9882, ...     2
    1.191e-3, ...   3
    0.8622, ...     4
    1.02501, ...    5
    1.0803e-4, ...  6
    0.99628, ...    7
    1.303e-4, ...   8
    1.02441, ...    9
    1.2986e-4, ...  10
    1.030, ...      11
    1.110e-2, ...   12
    1.10e-6, ...    13
    0.959, ...      14
    ];

% Tabata 1996a equations 13--19
c = [...
    D(1).*A ./ Z.^D(2), ... 1
    D(3).*Z.^D(4), ...      2
    D(5) - D(6).*Z, ...     3
    D(7) - D(8).*Z, ...     4
    D(9) - D(10).*Z, ...    5
    D(11)./Z.^D(12), ...    6
    D(13).*Z.^D(14), ...    7
    ];

% Tabata 1996a equation 12
B = log((t0 ./ (I + c(7).*t0)).^2) + log(1 + t0./2);

% Tabata 1996a equation 11
csdaRangeGCm2 = c(1) ./ B .* ...
    (log(1+c(2).*t0.^c(3))./c(2) - c(4).*t0.^c(5)./(1 + c(6).*t0));

