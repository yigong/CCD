% Points on the Plane of the Detector
N = 1000;  % Number of Points on Sphere
sr = 50; % Sphere Radius

        y0 = 3.0.*(rand(1,N)-0.5);
        x0 = 3.0.*(rand(1,N)-0.5);
%       // Point within one of the Detector's Layers
        z0 = - 0.650*(rand(1,N));

%       // Dirt Ball - Spherical Background Source
%       //G4double phi = 2*pi*G4UniformRand();
        cosTheta = (1-2.*(rand(1,N)));
        %// Semi-Isotropic Scattering
        %a =.5;
        % cosTheta = (2-1-4.*rand(1,N))/(1+sqrt(1-a.*(2-a-4.*rand(1,N))));
        %cosTheta2 = cosTheta
      
        phi = 2*pi*rand(1,N);  % ADDED an extra cos(theta) to account attenuation length increase at poles
        theta = acos(cosTheta);
        
        
        
        
%       // Cartiesain Points on 1m-radius sphere:
         x1 = sr.*cos(phi).*sin(theta);    
                 y1 = sr.*sin(phi).*sin(theta); % [cm]
                 z1 = sr.*cos(theta);           % [cm]

%        // Calculate Momentrum Direction
         x2 = (x0-x1); % [cm]
         y2 = (y0-y1); % [cm]
         z2 = (z0-z1); % [cm]

%        aMag = sqrt(x2*x2+y2*y2+z2*z2);
%        ux = (x2/aMag),
%                 uy = y2/aMag,
%                 uz = z2/aMag;
                 
                 plot3(x0,y0,z0,'.b', 'MarkerSize', 2); hold on
                 plot3(x2,y2,z2,'.r', 'MarkerSize',8)
%axis equal
xlabel('X');ylabel('Y');zlabel('Z')
axis equal