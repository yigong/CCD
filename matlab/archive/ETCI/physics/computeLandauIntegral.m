function omega = computeLandauIntegral(lambda_0)
% function omega = computeLandauIntegral(lambda_0)
% 
% Integrate Landau distribution above lambda_0
% as parameterized by D.H. Wilkinson, NIMA 383 (1996) 513-515. Accuracy is 0.01%

%% 1. Define parameters
% segment A: lambda_0 < -1
A0 = +0.008517;
A1 = +0.2892;
A2 = -0.8619;

% segment B: -1 <= lambda_0 <= 3
B0 = +0.713167;
B1 = -0.17885481;
B2 = +0.0077322340;
B3 = +0.010013627;
B4 = -0.0034453395;
B5 = +0.00022983112;
B6 = +0.00021393062;
B7 = -0.000089135946;
B8 = +0.000015782663;
B9 = -0.0000011232122;

% segment C: 3 <= lambda_0 <= 150
C0 = -0.50614708;
C1 = -0.59477519;
C2 = +0.33504975;
C3 = -0.35015962;
C4 = +0.083082286;
C5 = +0.01467629;
C6 = -0.011817844;
C7 = +0.0026055657;
C8 = -0.00026619565;
C9 = +0.000010770076;

% segment D: lambda_0 > 150
D0 = -4.1265648;
D1 = +3.6709034;
D2 = -1.1892577;
D3 = +0.13823403;
D4 = -0.0083406914;
D5 = +0.00020502525;
gammaE = 0.5772156649;  % euler constant

%% 2. Compute omega
% segment A: lambda_0 < -1
lgA = lambda_0 < -1;
lA = lambda_0(lgA);
omegaLess = 1 - (1 - erf(exp((abs(lA)-1)/2))) / sqrt(2);
omega(lgA) = omegaLess ./ (1 + A0*exp(A1*lA + A2*lA.^2));

% segment B: -1 <= lambda_0 <= 3
lgB = lambda_0 >= -1 & lambda_0 <= 3;
lB = lambda_0(lgB);
omega(lgB) = B0       + B1*lB    + B2*lB.^2 + B3*lB.^3 + B4*lB.^4 + ...
             B5*lB.^5 + B6*lB.^6 + B7*lB.^7 + B8*lB.^8 + B9*lB.^9;

% segment C: 3 <= lambda_0 <= 150
lgC = lambda_0 > 3 & lambda_0 <= 150;
LC = log(lambda_0(lgC));
omega(lgC) = exp(...
    C0       + C1*LC    + C2*LC.^2 + C3*LC.^3 + C4*LC.^4 + ...
    C5*LC.^5 + C6*LC.^6 + C7*LC.^7 + C8*LC.^8 + C9*LC.^9);

% segment D: lambda_0 > 150
lgD = lambda_0 > 150;
lD = lambda_0(lgD);
LD = log(lD);
omegaGreater = (lD + LD - 1 + gammaE) ./ lD.^2;
omega(lgD) = omegaGreater .* ...
    (1 + 0.01*exp(D0 + D1*LD + D2*LD.^2 + D3*LD.^3 + D4*LD.^4 + D5*LD.^5));

