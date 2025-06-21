function [A_lin, B_lin, C_lin, D_lin, K] = linearization()
    % Zdefiniowanie symboli
    syms theta dtheta r dr tau F real
    syms m g b1 b2 real

    % Stan i wejście
    x = [theta; dtheta; r; dr];
    u = [tau; F];

    % Równania dynamiki
    M = [m*r^2,         0;
         0,             m];

    C = [b1,            2*m*r*dtheta;
        -m*r*dtheta,    b2];

    G = [m*g*r*cos(theta);
         m*g*sin(theta)];

    v = [dtheta; dr];
    acc = M \ (u - C*v - G);

    dx = [dtheta; acc(1); dr; acc(2)];

    % Liniaryzacja symboliczna
    A_sym = jacobian(dx, x);
    B_sym = jacobian(dx, u);

    % Punkt pracy i parametry
    r0 = 1; 
    params = {theta,0; dtheta,0; r,r0; dr,0; m,1; g,9.81; b1,0.1; b2,0.1};

    % Zastąpienie symboli wartościami
    for i=1:length(params)
        A_sym = subs(A_sym, params{i,1}, params{i,2});
        B_sym = subs(B_sym, params{i,1}, params{i,2});
    end

    % Upraszczanie wyrażeń symbolicznych
    A_sym = simplify(A_sym);
    B_sym = simplify(B_sym);

    % Ręczne zastąpienie trygonometrii
    A_sym = subs(A_sym, {cos(0), sin(0)}, {1,0});
    B_sym = subs(B_sym, {cos(0), sin(0)}, {1,0});

    % Konwersja symboliczne → numeryczne
    A_lin = double(vpa(A_sym,8));
    B_lin = double(vpa(B_sym,8));

    % Sprawdzenie
    assert(isnumeric(A_lin) && isnumeric(B_lin), 'A lub B zawiera symbole!');

    % C i D (pełna obserwowalność)
    C_lin = eye(4);
    D_lin = zeros(4,2);

    % Regulator LQR
    Q = diag([100, 1, 100, 1]);
    R = eye(2);
    K = lqr(A_lin, B_lin, Q, R);

    % Wyświetlenie wyników
    disp('Macierz A:'); disp(A_lin);
    disp('Macierz B:'); disp(B_lin);
    disp('Macierz K (LQR):'); disp(K);
end
