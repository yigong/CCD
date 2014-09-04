% Points on the Plane of the Detector
N = 10000;  % Number of Points on Sphere
sr = 50; % Sphere Radius

for i=1:N
        y0(1,i) = 3.0.*(rand(1,1)-0.5);
        x0(1,i) = 3.0.*(rand(1,1)-0.5);
%       // Point within one of the Detector's Layers 
        z0(1,i) = - 0.650*(rand(1,1));

%       // Dirt Ball - Spherical Background Source
%       //G4double phi = 2*pi*G4UniformRand();
        % cosTheta(1,i) = (1-2.*(rand(1,1))); % For Isotropic Theta
        %theta(1,i) = acos(cosTheta(1,i));
        phi(1,i) = 2*pi*rand(1,1); 
        
        sinTheta(1,i) = sqrt(rand(1,1));
        theta(1,i) = asin(sinTheta(1,i));
        
        if rand(1) > 0.5
            theta(1,i) = pi - theta(1,i);
        end
end     
            
        
        
        
        
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