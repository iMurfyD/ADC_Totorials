clc; clear;
% gravitational parameter
G = 6.67e-11*(1e-9);
Me = 5.972e24; % Mass of the earth in kg
mu = G*Me; % Gravity

% Timespan
min = 60; hour = 60*min; day = 24*hour; year = 365*day;
t_start = 0;
dt = 1*min;
t_final = 3600*day;
timeSpan = t_start:dt:t_final;

% TLE data for ISS
i = 51.6440*(pi/180); % inclination in radians
e = .0006601; % eccentricity
M0 = 96.5029*(pi/180); % Mean anomally in rad
Omega_RAAN = 214.1057*(pi/180); % Right asscention in radians
freq = 15.54082521; % rev's per day
w = 54.8109*(pi/180); % Argument of periapsis
T = (freq*(1/24)*(1/60)*(1/60))^-1; % Period of orbit (s) 


a = nthroot((mu*(T^2))/(4*pi^2),3); % semi-major axis (km)

[r0, v0] = elements_to_sv(a,e,i,Omega_RAAN,w,M0,dt); 
X0 = [r0; v0];
 
[t_out, X_out] = ode45(@twobody,timeSpan,X0, odeset(),mu);

r = zeros(3,length(t_out));

for t = 1:length(t_out)
    r(1:3,t) = X_out(t,1:3);

end

PlotOrbits(r);


