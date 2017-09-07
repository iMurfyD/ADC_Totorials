clc; clear all;

% gravitational parameter
G = 6.67e-11*(1e-9);
Me = 5.972e24; % Mass of the earth in kg
mu = G*Me; % Gravity


% Set timespan
sec = 1;min = 60*sec; hours = 60*min;
ti = 0*sec;
dt = 10*sec; % Timestep of 10s
tf = 10*min;
timeSpan = ti:dt:tf; % Array of timesteps

% Dimensions of our .5x.5x.5 m box
w = .5; % width = 1/2 m
d = .5; % depth = 1/2 m
h = .5; % height = 1/2 m
m = 12; % 12 kg mass
J = [(1/12)*m*(w^2+h^2) 0 0; % Inertia matrix in kg*m^2
     0 (1/12)*m*(d^2+w^2) 0;
     0 0 (1/12)*m*(d^2+h^2);];
 
ang_w = [1 0 0]'*(pi/180); % in deg -> rad/s

q0 = [ .999 .002 .628 -.209]'; % From quaternions.online
[t_out, q_out] = ode45(@qdot,timeSpan,q0, odeset(), ang_w);

A_matrix = zeros(length(t_out),9);
% Fill in attitude matrix with quaternions
for t = 1:length(t_out)
    q = q_out(t,1:end);
    A_matrix(t,1) = q(1)^2-q(2)^2-q(3)^2+q(4)^2; % a(1,1)
    A_matrix(t,2) = 2*(q(1)*q(2)+q(3)*q(4)); % a(1,2)
    A_matrix(t,3) = 2*(q(1)*q(3)-q(2)-q(4)); % a(1,3)
    A_matrix(t,4) = 2*(q(1)*q(2)-q(3)*q(4)); % a(2,1)
    A_matrix(t,5) = -q(1)^2+q(2)^2-q(3)^2+q(4)^2; % a(2,2)
    A_matrix(t,6) = 2*(q(2)*q(3)+q(1)*q(4)); % a(2,3)
    A_matrix(t,7) = 2*(q(1)*q(3)-q(2)*q(4)); % a(3,1)
    A_matrix(t,8) = 2*(q(2)*q(3)-q(1)*q(4)); % a(3,2)
    A_matrix(t,9) = -q(1)^2-q(2)^2+q(3)^2+q(4)^2; % a(3,3)
end

% Make figure larger
hFig = figure(1);
set(hFig, 'Position', [500 500 1000*(1.3) 800*(1.3)]);

% Plot stuff
plot( t_out,q_out(1:end,1),'r', ...
      t_out,q_out(1:end,2),'g', ...
      t_out,q_out(1:end,3),'k', ...
      t_out,q_out(1:end,4),'b');
xlabel('Time (s)');
ylabel('Quaternion');
legend('q1','q2','q3','q4');

