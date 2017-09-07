function dAdt = angular_vel_inv_J(t, A, J, L)
dAdt = zeros(6,1);
dAdt(1:3) = A(4:6);
dAdt(4:6) = inv(J)*L;
end