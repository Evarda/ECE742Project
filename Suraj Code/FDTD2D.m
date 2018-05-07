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
ep(1:maxLength,  1:maxLength)   = 0;
mu(1:maxLength-1,1:maxLength)   = 0;

boundSize=10; %also called d?
% ep(1,boundSize,:)=boundary decay value
% ep(maxlangth-boundasize,maxlength,:)=boundary decay value

% S-matrix
% kx=
% ky=
% kz=1;
% sigmax=
% sigmay=
% sigmaz=0;
% sx=kx+sigmax/(i*2*pi*f*ep0);
% sy=ky+sigmay/(i*2*pi*f*ep0);
% sz=kz+sigmaz/(i*2*pi*f*ep0);

% Constants to update D
% CDx1=(2*ep0*ky-sigmay*dt)/(2*ep0*ky+sigmay*dt);
% CDx2=(2*ep0*dt)/(2*ep0*ky+sigmay*dt);
% CDy1=(2*ep0*kz-sigmaz*dt)/(2*ep0*kz+sigmaz*dt);
% CDy2=(2*ep0*dt)/(2*ep0*kz+sigmaz*dt);
% CDz1=(2*ep0*kx-sigmax*dt)/(2*ep0*kx+sigmax*dt);
% CDz2=(2*ep0*dt)/(2*ep0*kx+sigmax*dt);

% Constants to update E 
% CEx1=(2*ep0*kz-sigmaz*dt)/(2*ep0*kz+sigmaz*dt);
% CEx2=1/(2*ep0*kz+sigmaz*dt);
% CEy1=(2*ep0*kx-sigmax*dt)/(2*ep0*kx+sigmax*dt);
% CEy2=1/(2*ep0*kx+sigmax*dt);
% CEz1=(2*ep0*ky-sigmay*dt)/(2*ep0*ky+sigmay*dt);
% CEz2=1/(2*ep0*ky89+sigmay*dt);

% Constants to update E (used for D)
% CEDx1=2*ep0*kx+sigmax*dt;
% CEDx2=2*ep0*kx-sigmax*dt;
% CEDy1=2*ep0*ky+sigmay*dt;
% CEDy2=2*ep0*ky-sigmay*dt;
% CEDz1=2*ep0*kz+sigmaz*dt;
% CEDz2=2*ep0*kz-sigmaz*dt;

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

% Update loop for UPML region
