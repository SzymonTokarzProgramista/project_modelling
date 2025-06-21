function dx = nonlinear_dynamics(x, u, m, g, b1, b2)
% Oblicza pochodną stanu dla ramienia w układzie biegunowym
% x = [theta; r; dtheta; dr]
% u = [tau; F]

    theta = x(1);
    r = x(2);
    dtheta = x(3);
    dr = x(4);
    tau = u(1);
    F = u(2);

    % Równania ruchu:
    ddtheta = (tau - 2 * m * r * dr * dtheta - b1 * dtheta) / (m * r^2);
    ddr = (F + m * r * dtheta^2 - m * g * cos(theta) - b2 * dr) / m;

    dx = [dtheta; dr; ddtheta; ddr];
end
