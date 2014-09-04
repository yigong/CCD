function phi = computeLandau(lambda)
% function phi = computeLandau(lambda)
% 
% Sample Landau distribution 
% as parameterized by D.H. Wilkinson, NIMA 383 (1996) 513-515. Accuracy is 0.01%

% % e: particle charge
% % v: particle velocity
% % rho: density of absorber
% % Z, A: atomic number and weight of absorber
% Nav = 6.022e23; % Avogadro's number
% m = 511;    % mass of electron
% % x = thickness of absorber traversed
% 
% omega = (2*pi * Nav *e^4 * rho / m / v^2) * (Z / A) / (eps^2);
% 
% % beta = velocity / c
% % S = x / rho [mg/cm^2]
% eta = (0.1536 / beta^2) * (Z/A) * S;
% % eta also defined exactly using factors from omega

% segmented parameterization
%% 1. define parameters
%   segment A: -1 < lambda < 3
A0 = +0.17885481;
A1 = -0.015464468;
A2 = -0.030040882;
A3 = +0.013781358;
A4 = -0.0011491556;
A5 = -0.0012835837;
A6 = +0.00062395162;
A7 = -0.0001262613;
A8 = +0.000010108918;

%   segment B: 3 < lambda < 150
B0 = -1.5669876;
B1 = -1.5811676;
B2 = +1.677088;
B3 = -1.4972908;
B4 = +0.57062974;
B5 = -0.11777036;
B6 = +0.01343737;
B7 = -0.00076315158;
B8 = +0.000014881108;

%   segment C: lambda > 150
C0 = 5.157;
C1 = -1.42;
gammaE = 0.5772156649;  % euler constant

%   segment D: -3.4 < lambda < -1
D0 = 6.7853;
D1 = 4.884;
D2 = 1.4488;
D3 = 0.20802;
D4 = 0.012057;

%   segment E: lambda < -3.4 (inaccurate and unimportant)

%% 2. Calculate phi
%   segment A
lgA = lambda >= -1 & lambda <= 3;
lA = lambda(lgA);
phi(lgA) = A0 + A1*lA    + A2*lA.^2 + A3*lA.^3 + A4*lA.^4 + ...
                A5*lA.^5 + A6*lA.^6 + A7*lA.^7 + A8*lA.^8;
%   segment B
lgB = lambda > 3 & lambda <= 150;
LB = log(lambda(lgB));
phi(lgB) = exp(...
    B0 + B1*LB    + B2*LB.^2 + B3*LB.^3 + B4*LB.^4 + ...
         B5*LB.^5 + B6*LB.^6 + B7*LB.^7 + B8*LB.^8);

%   segment C
lgC = lambda > 150;
lC = lambda(lgC);
LC = log(lC);
phiGreater = lC.^-2 - (3 - 2*gammaE - 2*LC).*lC.^-3;
phi(lgC) = phiGreater ./ (1 - 0.01*exp(C0 + C1*LC));

%   segment D
lgD = lambda >= -3.4 & lambda < -1;
lD = lambda(lgD);
phiLess = exp( (abs(lD)-1)/2 - exp(abs(lD)-1) ) / sqrt(2*pi);
phi(lgD) = phiLess .* (1 + 0.01*(D0 + D1*lD + D2*lD.^2 + D3*lD.^3 + D4*lD.^4));

%   segment E
lgE = lambda < -3.4;
if any(lgE)
    warning('Less accurate for lambda < -3.4, but this region is insignificant.')
    lE = lambda(lgE);
    phiLess = exp( (abs(lE)-1)/2 - exp(abs(lE)-1)) / sqrt(2*pi);
    phi(lgE) = phiLess;
end
