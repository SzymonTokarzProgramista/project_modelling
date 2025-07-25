%parameters
g=9.8036; %m/s^2
m=0.01; %kg
rmin=0.1; %m
rmax=0.4; %m
theta_min = -pi/2;%rad
theta_max=pi/2;%rad
dtheta_min=-10;%rad/s
dtheta_max=10;%rad/s
dr_min=-10;%m/s
dr_max=10;%m/s
b1=0.2;
b2=0.2;


%initial conditions
theta_start = pi/6;%rad
r_start = 0.125; %m


% Parametry trajektorii
theta_start = pi/6;
theta_desired = -pi/3;

r_start = 0.125;
r_desired = 0.3;

t = [0 5 10 15 20]; % czas zmiany pozycji [sekundy]
theta_ref = [theta_start  theta_desired  theta_desired  theta_start theta_start]; % kąt najpierw 0, potem pi/4, potem powrót do 0
r_ref = [r_start r_desired r_desired r_start r_start];    % długość ramienia: najpierw krótka, potem długa, powrót

% Dane do bloku From Workspace:
theta_input = [t' theta_ref'];
r_input = [t' r_ref'];

% Liczba punktów
N = 20;

% Interpolacja punktów w przestrzeni biegunowej
theta_points = linspace(theta_start, theta_desired, N);
r_points = linspace(r_start, r_desired, N);

% Inicjalizacja wyników
pid_theta = zeros(N, 3); % [Kp Ki Kd]
pid_r = zeros(N, 3);
xy_points = zeros(N, 2);

for i = 1:N
    theta = theta_points(i);
    r = r_points(i);

    % Zamiana na współrzędne robocze
    x = r * cos(theta);
    y = r * sin(theta);
    xy_points(i,:) = [x, y];

    % Wywołanie funkcji strojącej PIDy
    pid = compute_pid_from_target(x, y);

    pid_theta(i,:) = [pid.theta.Kp, pid.theta.Ki, pid.theta.Kd];
    pid_r(i,:) = [pid.r.Kp, pid.r.Ki, pid.r.Kd];
end

% Wyświetlenie
disp('PIDy dla trajektorii:');
T = table(xy_points(:,1), xy_points(:,2), ...
    pid_theta(:,1), pid_theta(:,2), pid_theta(:,3), ...
    pid_r(:,1), pid_r(:,2), pid_r(:,3), ...
    'VariableNames', {'x','y','Kp_theta','Ki_theta','Kd_theta','Kp_r','Ki_r','Kd_r'});
disp(T);


theta_grid = linspace(-pi/2, pi/2, 20);
r_grid = linspace(0.1, 0.4, 20);

[Kp_theta_map, Ki_theta_map, Kd_theta_map] = deal(zeros(20,20));
[Kp_r_map, Ki_r_map, Kd_r_map] = deal(zeros(20,20));

for i = 1:20
    for j = 1:20
        theta = theta_grid(i);
        r = r_grid(j);

        x = r * cos(theta);
        y = r * sin(theta);

        pid = compute_pid_from_target(x, y);

        % Ograniczenie theta-PID
        pid.theta.Kp = min(max(pid.theta.Kp, 0), 100);     % [0, 100]
        pid.theta.Ki = min(max(pid.theta.Ki, 0), 10);      % [0, 10]
        pid.theta.Kd = min(max(pid.theta.Kd, 0), 5);       % [0, 5]

        pid.r.Kp = min(max(pid.r.Kp, 0), 100);     % [0, 100]
        pid.r.Ki = min(max(pid.r.Ki, 0), 10);      % [0, 10]
        pid.r.Kd = min(max(pid.r.Kd, 0), 5);       % [0, 5]



        Kp_theta_map(i,j) = pid.theta.Kp;
        Ki_theta_map(i,j) = pid.theta.Ki;
        Kd_theta_map(i,j) = pid.theta.Kd;

        Kp_r_map(i,j) = pid.r.Kp;
        Ki_r_map(i,j) = pid.r.Ki;
        Kd_r_map(i,j) = pid.r.Kd;
    end
end

save('pid_maps.mat', ...
    'Kp_theta_map','Ki_theta_map','Kd_theta_map', ...
    'Kp_r_map','Ki_r_map','Kd_r_map', ...
    'theta_grid','r_grid');
