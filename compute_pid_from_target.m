function pid_gains = compute_pid_from_target(x_target, y_target)
    % Parametry układu
    g = 9.8036;
    m = 0.01;
    b1 = 0.2;
    b2 = 0.2;

    % 1. Przekształcenie (x,y) -> (r, theta)
    r = hypot(x_target, y_target);
    theta = atan2(y_target, x_target);

    % Stan początkowy
    dtheta = 0;
    dr = 0;
    x = [theta; r; dtheta; dr];
    u = [0; 0]; % brak wejścia do liniaryzacji

    % Nieliniowa dynamika
    dyn = @(x,u) nonlinear_dynamics(x,u,m,g,b1,b2);

    % 2. Liniaryzacja numeryczna
    delta = 1e-6;
    A = zeros(4);
    B = zeros(4,2);

    for i = 1:4
        dx = zeros(4,1); dx(i) = delta;
        A(:,i) = (dyn(x + dx, u) - dyn(x - dx, u)) / (2*delta);
    end
    for i = 1:2
        du = zeros(2,1); du(i) = delta;
        B(:,i) = (dyn(x, u + du) - dyn(x, u - du)) / (2*delta);
    end

    % 3. Wyodrębnienie podukładów
    A_theta = A([1 3],[1 3]);
    B_tau   = B([1 3],1);

    A_r = A([2 4],[2 4]);
    B_F = B([2 4],2);

    % 4. Tworzenie modeli SS i strojenie PID
    sys_theta = ss(A_theta, B_tau, [1 0], 0);
    sys_r = ss(A_r, B_F, [1 0], 0);

    C_theta = pidtune(sys_theta, 'PID');
    C_r = pidtune(sys_r, 'PID');

    % 5. Wynik
    pid_gains = struct();
    pid_gains.theta = struct('Kp', C_theta.Kp, 'Ki', C_theta.Ki, 'Kd', C_theta.Kd);
    pid_gains.r     = struct('Kp', C_r.Kp,     'Ki', C_r.Ki,     'Kd', C_r.Kd);

    % Wypisz
    disp('PID dla τ (sterowanie kątem θ):');
    disp(pid_gains.theta);
    disp('PID dla F (sterowanie promieniem r):');
    disp(pid_gains.r);
end
