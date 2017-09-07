function dXdt = twobody(t, X, mu) % X = [r v]

dXdt = zeros(6, 1);
dXdt(1:3) = X(4:6);
dXdt(4:6) = -mu .* X(1:3) ./ (norm(X(1:3)).^3);
%fprintf('%g\n',t);
end
   
