% Material Grid Properties
% Uniform Free Space
mu0 = 4*pi*10^-7; %[H/m]
ep0 = 8.854187817*10^-12; %[F/m]

c = 299792458; % speed of light [m/s]
f = 10^9; % frequency [1/s]
lambda = c/f; % wavelength [m]
dx = lambda/10; % step size x [m]
dy = lambda/10; % step size y [m]
dt = dx/(sqrt(2)*c); % step size t

% Declare Grids
maxLength=50;
[iarray,jarray] = meshgrid(1:maxLength,1:maxLength);
Ez(1:maxLength,  1:maxLength)   = 0;
Hx(1:maxLength,  1:maxLength-1) = 0;
Hy(1:maxLength-1,1:maxLength)   = 0;

% Calculate Update Factor
HxUpFacy = dt/(mu0*dy);
HyUpFacx = dt/(mu0*dx);
EzUpFacx = dt/(ep0*dx);
EzUpFacy = dt/(ep0*dy);

% Max Timesteps
nmax = 200;

% Initialize figure
figure

% Update Loop
for n = 1:nmax
    for i = 1:maxLength
        for j=1:maxLength-1
        % Update Hx
        Hx(i,j)=Hx(i,j)-HxUpFacy*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1:maxLength-1
        for j=1:maxLength
        % Update Hy
        Hy(i,j)=Hy(i,j)+HyUpFacx*(Ez(i+1,j)-Ez(i,j));
        end
    end
    for i = 2:maxLength-1
        for j=2:maxLength-1
        % Update Ez
        Ez(i,j)=Ez(i,j)+(EzUpFacx*(Hy(i,j)-Hy(i-1,j))-EzUpFacy*(Hx(i,j)-Hx(i,j-1)));
        end
    end
    
    % PEC Boundary Conditions
    Ez(1:maxLength,1) = 0;
    Ez(1:maxLength, maxLength) = 0;
    Ez(1,1:maxLength) = 0;
    Ez(maxLength, 1:maxLength) = 0;
    
    % Source
    Ez(maxLength/2,maxLength/2) = cos(2*pi*f*dt*n);
    
    % Plot
    % Surface Plot
    subplot(2,1,1)
    surf(iarray, jarray, Ez);
    axis([1, maxLength, 1, maxLength, -1, 1])
    xlabel('x [m]')
    xticks(iarray(1,10:10:maxLength))
    xticklabels(dx*iarray(1,10:10:maxLength))
    ylabel('y [m]')
    yticks(jarray(10:10:maxLength,1))
    yticklabels(dy*jarray(10:10:maxLength,1))
    zlabel('E')
    
    % Contour Plot
    subplot(2,1,2)
    contour(iarray, jarray, Ez);
    xlabel('x [m]')
    xticks(iarray(1,10:10:maxLength))
    xticklabels(dx*iarray(1,10:10:maxLength))
    ylabel('y [m]')
    yticks(jarray(10:10:maxLength,1))
    yticklabels(dy*jarray(10:10:maxLength,1))
    zlabel('E')
    M=getframe;
    
end