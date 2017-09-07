% Converts orbital elements to a state vector.
% By Benjamin Reifler, last updated 2016-04-09.

function [r,v]=elements_to_sv(a,e,i,Omega,w,M0,dt)
mu = 398600.64;
n = sqrt(mu/a^3);
M = M0 + n*dt;

while M >= 2*pi
    M = M-2*pi;
end
while M < 2*pi
    M = M+2*pi;
end

% Use Newton-Raphson method to solve for E
% f(E) : 0 = E-e*sin(E)-M
f = @(E) E-e*sin(E)-M;
% df(E)/dE : 0 = 1-e*cos(E)
df = @(E) 1-e*cos(E);

E = 1;
%n = 0;
F = 1;
c = 1;
eps = 0.001;
while abs(F) > eps && abs(c) > eps
    %n = n + 1;
    F = f(E);
    dF = df(E);
    c = -F/dF;
    E = E + c;
end

rmag = a*(1-e*cos(E));

x = a*(cos(E)-e);
y = a*sqrt(1-e^2)*sin(E);
xdot = -sqrt(mu*a)/rmag*sin(E);
ydot = sqrt(mu*a*(1-e^2))/rmag*cos(E);

% Attitude matrix ECI -> perifocal
so = sin(Omega);
co = cos(Omega);
sw = sin(w);
cw = cos(w);
si = sin(i);
ci = cos(i);

A = [co*cw-so*sw*ci  so*cw+co*sw*ci  sw*si;...
     -co*sw-so*cw*ci -so*sw+co*cw*ci cw*si;...
     so*si           -co*si          ci];

r = A'*[x y 0]';
v = A'*[xdot ydot 0]';