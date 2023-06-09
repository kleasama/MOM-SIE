function [ k,E,H ] = PW( amplitude,wavenumber,eta,theta,phi,alpha )
%PlaneWave Create a time-harmonic electromagnetic plane wave.

%% k components
kx = sin(theta)*cos(phi);
ky = sin(theta)*sin(phi);
kz = cos(theta);
%% E components
Ex = cos(theta)*cos(phi)*cos(alpha) - sin(phi)*sin(alpha);
Ey = cos(theta)*sin(phi)*cos(alpha) + cos(phi)*sin(alpha);
Ez = -sin(theta)*cos(alpha);
%% H components
Hx = cos(theta)*cos(phi)*cos(alpha + pi/2) - sin(phi)*sin(alpha + pi/2);
Hy = cos(theta)*sin(phi)*cos(alpha + pi/2) + cos(phi)*sin(alpha + pi/2);
Hz = -sin(theta)*cos(alpha + pi/2);
%% k, E, H vectors
k  = wavenumber*[kx,ky,kz];
E  = amplitude*[Ex,Ey,Ez];
H  = amplitude/eta*[Hx,Hy,Hz];
end