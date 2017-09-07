function dqdt = qdot(t, q, w)
% Taken from Crassidis larc lecture on attitude determination and
% estimation
% https://mediaex-server.larc.nasa.gov/Academy/Play/bdeb764e048940a6b2ae05c3cfdf5d261d?catalog=8e500782-c73d-4bc2-ad69-cf59aec8420c

wcross = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
Omega = [-wcross w; -w' 0];
dqdt = .5*Omega*q;
end