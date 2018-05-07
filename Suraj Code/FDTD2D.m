clc;
close all;
clear all;

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
maxLength=100;
[iarray,jarray] = meshgrid(1:maxLength,1:maxLength);
Ez(1:maxLength,  1:maxLength)   = 0;
Dz(1:maxLength,  1:maxLength)   = 0;
Hx(1:maxLength,  1:maxLength) = 0;
Bx(1:maxLength,  1:maxLength) = 0;
Hy(1:maxLength, 1:maxLength)   = 0;
By(1:maxLength, 1:maxLength)   = 0;

ep=1;
mu=1;
boundsize=10; %also called d?
% ep(1,boundSize,:)=boundary decay value
% ep(maxlangth-boundasize,maxlength,:)=boundary decay value

% S-matrix
kx=1;
ky=1;
kz=1;
sigmax=0;
sigmay=0;
sigmaz=0;
% sigmax(1:boundsize)=0;
% sigmay(1:boundsize)=0;
% sigmaz(1:boundsize)=0;
% sx=kx+sigmax/(i*2*pi*f*ep0);
% sy=ky+sigmay/(i*2*pi*f*ep0);
% sz=kz+sigmaz/(i*2*pi*f*ep0);

% Constants to update B
CBx1=(2*ep0*ky-sigmay*dt)/(2*ep0*ky+sigmay*dt);
CBx2=(2*ep0*dt)/(2*ep0*ky+sigmay*dt);
CBy1=(2*ep0*kz-sigmaz*dt)/(2*ep0*kz+sigmaz*dt);
CBy2=(2*ep0*dt)/(2*ep0*kz+sigmaz*dt);
% Constants to update D
CDz1=(2*ep0*kx-sigmax*dt)/(2*ep0*kx+sigmax*dt);
CDz2=(2*ep0*dt)/(2*ep0*kx+sigmax*dt);

% Constants to update H
CHx1=(2*ep0*kz-sigmaz*dt)/(2*ep0*kz+sigmaz*dt);
CHx2=1/(2*ep0*kz+sigmaz*dt);
CHy1=(2*ep0*kx-sigmax*dt)/(2*ep0*kx+sigmax*dt);
CHy2=1/(2*ep0*kx+sigmax*dt);
% Constants to update E
CEz1=(2*ep0*ky-sigmay*dt)/(2*ep0*ky+sigmay*dt);
CEz2=1/(2*ep0*ky+sigmay*dt);

% Constants to update H (used for B)
CHBx1=2*ep0*kx+sigmax*dt;
CHBx2=2*ep0*kx-sigmax*dt;
CHBy1=2*ep0*ky+sigmay*dt;
CHBy2=2*ep0*ky-sigmay*dt;
% Constants to update E (used for D)
CEDz1=2*ep0*kz+sigmaz*dt;
CEDz2=2*ep0*kz-sigmaz*dt;

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
    % Update Hx
    for i = 1:boundsize-1
        for j=1:boundsize-1
            Bx_old=Bx(i,j);
            Bx(i,j)=CBx1*Bx(i,j)+CBx2*(Ez(i,j+1)-Ez(i,j))/dy;
            Hx(i,j)=CHx1*Hx(i,j)+CHx2*(CHBx1*Bx(i,j)-CHBx2*Bx_old)/mu;
        end
    end    
    for i = boundsize:maxLength-boundsize-1
        for j=boundsize:maxLength-boundsize-1
            Hx(i,j)=Hx(i,j)-HxUpFacy*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = maxLength-boundsize:maxLength-1
        for j=maxLength-boundsize:maxLength-1
            Bx_old=Bx(i,j);
            Bx(i,j)=CBx1*Bx(i,j)+CBx2*(Ez(i,j+1)-Ez(i,j))/dy;
            Hx(i,j)=CHx1*Hx(i,j)+CHx2*(CHBx1*Bx(i,j)-CHBx2*Bx_old)/mu;
        end
    end
    
    % Update Hy
    for i = 1:boundsize-1
        for j=1:boundsize-1
            By_old=By(i,j);
            By(i,j)=CBy1*Bx(i,j)-CBy2*(Ez(i,j+1)-Ez(i,j))/dx;
            Hy(i,j)=CHy1*Hx(i,j)+CHy2*(CHBy1*Bx(i,j)-CHBy2*By_old)/mu;
        end
    end    
    for i = boundsize:maxLength-boundsize-1
        for j=boundsize:maxLength-boundsize-1         
            Hy(i,j)=Hy(i,j)+HyUpFacx*(Ez(i+1,j)-Ez(i,j));
        end
    end
    for i = maxLength-boundsize:maxLength-1
        for j=maxLength-boundsize:maxLength-1
            By_old=By(i,j);
            By(i,j)=CBy1*Bx(i,j)-CBy2*(Ez(i,j+1)-Ez(i,j))/dx;
            Hy(i,j)=CHy1*Hx(i,j)+CHy2*(CHBy1*Bx(i,j)-CHBy2*By_old)/mu;
        end
    end
    
    % Update Ez  
    for i = 2:boundsize-1
        for j=2:boundsize-1
            Dz_old=Dz(i,j);
            Dz(i,j)=CDz1*Dz(i,j)+CDz2*((Hy(i,j)-Hy(i-1,j))/dx-(Hx(i,j)-Hx(i,j-1))/dy);
            Ez(i,j)=CEz1*Ez(i,j)+CEz2*(CEDz1*Dz(i,j)-CEDz2*Dz_old)/ep;
        end
    end
    for i = boundsize:maxLength-boundsize-1
        for j=boundsize:maxLength-boundsize-1
            Ez(i,j)=Ez(i,j)+(EzUpFacx*(Hy(i,j)-Hy(i-1,j))-EzUpFacy*(Hx(i,j)-Hx(i,j-1)));
        end
    end
    for i = maxLength-boundsize:maxLength
        for j=maxLength-boundsize:maxLength
            Dz_old=Dz(i,j);
            Dz(i,j)=CDz1*Dz(i,j)+CDz2*((Hy(i,j)-Hy(i-1,j))/dx-(Hx(i,j)-Hx(i,j-1))/dy);
            Ez(i,j)=CEz1*Ez(i,j)+CEz2*(CEDz1*Dz(i,j)-CEDz2*Dz_old)/ep;
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
