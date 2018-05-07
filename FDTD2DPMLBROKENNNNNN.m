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
dt = dx/(c*sqrt(2)); % step size t

% Declare Grids
maxLength=100;
[iarray,jarray] = meshgrid(1:maxLength,1:maxLength);
Ez(1:maxLength,  1:maxLength)   = 0;
Dz(1:maxLength,  1:maxLength)   = 0;
Hx(1:maxLength,  1:maxLength) = 0;
Bx(1:maxLength,  1:maxLength) = 0;
Hy(1:maxLength, 1:maxLength)   = 0;
By(1:maxLength, 1:maxLength)   = 0;

ep=1*ep0;
mu=1*mu0;
boundsize=20;
% ep(1,boundSize,:)=boundary decay value
% ep(maxlangth-boundasize,maxlength,:)=boundary decay value

% S-matrix Values
kx(1:maxLength)=1;
ky(1:maxLength)=1;
kz(1:maxLength)=1;
sigmax(1:maxLength)=0;
sigmay(1:maxLength)=0;
sigmaz(1:maxLength)=0;

% Polynomial Grading
m=2;
sigmamax=5;
for i=1:boundsize
    sigmax(i)=((i/boundsize)^m*sigmamax);
    sigmay(i)=((i/boundsize)^m*sigmamax);
    sigmax(maxLength-boundsize+i)=sigmax(i);
    sigmay(maxLength-boundsize+i)=sigmay(i);
end

% Declare Arrays

% Constants to Update Hx
CBX1(1:maxLength)=0;
CBX2(1:maxLength)=0;
CHX1(1:maxLength)=0;
CHX2(1:maxLength)=0;
CHX3(1:maxLength)=0;
CHX4(1:maxLength)=0;

% Constants to Update Hy
CBY1(1:maxLength)=0;
CBY2(1:maxLength)=0;
CHY1(1:maxLength)=0;
CHY2(1:maxLength)=0;
CHY3(1:maxLength)=0;
CHY4(1:maxLength)=0;

% Constants to update Ez
CDZ1(1:maxLength)=0;
CDZ2(1:maxLength)=0;
CEZ1(1:maxLength)=0;
CEZ2(1:maxLength)=0;
CEZ3(1:maxLength)=0;
CEZ4(1:maxLength)=0;



% Fill Arrays

for i=1:maxLength
    
% Constants to Update Hx - index j j 0 0 i i
CBX1(i)=(2*ep0*ky(i)-sigmay(i)*dt)/(2*ep0*ky(i)+sigmay(i)*dt);
CBX2(i)=(2*ep0*dt)/(2*ep0*ky(i)+sigmay(i)*dt);
CHX1(i)=(2*ep0*kz(i)-sigmaz(i)*dt)/(2*ep0*kz(i)+sigmaz(i)*dt); %Constant
CHX2(i)=1/(2*ep0*kz(i)+sigmaz(i)*dt); %Constant
CHX3(i)=2*ep0*kx(i)+sigmax(i)*dt;
CHX4(i)=2*ep0*kx(i)-sigmax(i)*dt;

% Constants to Update Hy - index 0 0 i i j j
CBY1(i)=(2*ep0*kz(i)-sigmaz(i)*dt)/(2*ep0*kz(i)+sigmaz(i)*dt);
CBY2(i)=(2*ep0*dt)/(2*ep0*kz(i)+sigmaz(i)*dt); %Constant
CHY1(i)=(2*ep0*kx(i)-sigmax(i)*dt)/(2*ep0*kx(i)+sigmax(i)*dt);
CHY2(i)=1/(2*ep0*kx(i)+sigmax(i)*dt);
CHY3(i)=2*ep0*ky(i)+sigmay(i)*dt;
CHY4(i)=2*ep0*ky(i)-sigmay(i)*dt;

% Constants to update Ez - index i i j j 0 0
CDZ1(i)=(2*ep0*kx(i)-sigmax(i)*dt)/(2*ep0*kx(i)+sigmax(i)*dt);
CDZ2(i)=(2*ep0*dt)/(2*ep0*kx(i)+sigmax(i)*dt);
CEZ1(i)=(2*ep0*ky(i)-sigmay(i)*dt)/(2*ep0*ky(i)+sigmay(i)*dt);
CEZ2(i)=1/(2*ep0*ky(i)+sigmay(i)*dt);
CEZ3(i)=2*ep0*kz(i)+sigmaz(i)*dt;
CEZ4(i)=2*ep0*kz(i)-sigmaz(i)*dt;

end

% Max Timesteps
nmax = 200;

% Initialize figure
figure

% Update Loop
for n = 1:nmax
    % Update Hx
    for i = 1:maxLength-1
        for j=1:maxLength-1
            Bx_old=Bx(i,j);
            Bx(i,j)=CBX1(j)*Bx(i,j)+CBX2(j)*(Ez(i,j+1)-Ez(i,j))/dy;
            Hx(i,j)=CHX1(0)*Hx(i,j)+CHX2(0)*(CHX3(i)*Bx(i,j)-CHX4(i)*Bx_old)/mu;
        end
    end    

    
    % Update Hy
    for i = 1:maxLength-1
        for j=1:maxLength-1
            By_old=By(i,j);
            By(i,j)=CBY1(0)*Bx(i,j)-CBY2(0)*(Ez(i,j+1)-Ez(i,j))/dx;
            Hy(i,j)=CHY1(i)*Hx(i,j)+CHY2(i)*(CHY3(j)*Bx(i,j)-CHY4(j)*By_old)/mu;
        end
    end    
    
    % Update Ez  
    for i = 2:maxLength
        for j=2:maxLength
            Dz_old=Dz(i,j);
            Dz(i,j)=CDZ1(i)*Dz(i,j)+CDZ2(i)*((Hy(i,j)-Hy(i-1,j))/dx-(Hx(i,j)-Hx(i,j-1))/dy);
            Ez(i,j)=CEZ1(j)*Ez(i,j)+CEZ2(j)*(CEZ3(0)*Dz(i,j)-CEZ4(0)*Dz_old)/ep;
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
