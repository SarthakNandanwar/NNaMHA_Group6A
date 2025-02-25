close all
clear
%%Single Cell
%define parameters
g_0 = 5e-3;
g_1 = 1700;
C_th = 1e-14;
g_th= 1/1.5e6;
T_0 = 300;
V_b = 1.1;
R_b = 1e3;
C = 10e-9;

%Initial conditions
T_init = 300;
V_c_init = 0.001;


t_span = [0 : 0.001e-4 : 2.5e-4];

G = @(T) g_0 .* exp(-g_1 ./ T)

%Define system of ODEs
ode_system = @(t, y) [((V_b - y(1)) ./ R_b - G(y(2)) .* y(1)) ./ C; 
                      ((y(1) .^ 2) .* G(y(2)) - (y(2)-T_0) .* g_th) ./ C_th]
% Initial conditions vector [V_c, T]
initial_conditions = [V_c_init; T_init];

% Solve the system using ode45
[t, y] = ode45(ode_system, t_span, initial_conditions);

% Extract the solutions for V_c and T
V_c = y(:, 1); % Voltage across the capacitor
T = y(:, 2);   % Temperature



figure;
subplot(1,2,1)
plot(t, V_c, 'b-', 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Voltage across capacitor V_c(t)');
grid on;

subplot(1,2,2)
plot(t, T, 'r-', 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Temperature T(t)');
title('Temperature vs. Time');
grid on;


%%Single Cell with Rs
%define parameters
g_0 = 5e-3;
g_1 = 1700;
C_th = 1e-14;
g_th= 1/1.5e6;
T_0 = 300;
V_b = 1.1;
R_b = 1e3;
C = 10e-9;
R_s = 200;

%Initial conditions
T_init = 300;
V_c_init = 0.001;

%Timespan
t_span2 = [0 : 0.00001e-4 : 2.5e-4];

%Functions
G = @(T) g_0 .* exp(-g_1 ./ T);
V_mc = @(V_c, T) V_c ./ (1 + G(T) .* R_s);
%Define system of ODEs
ode_system = @(t, y) [((V_b - y(1)) ./ R_b - G(y(2)) .* y(1)) ./ C; 
                      ((V_mc(y(1), y(2)) .^ 2) .* G(y(2)) - (y(2)-T_0) .* g_th) ./ C_th]
% Initial conditions vector [V_c, T]
initial_conditions = [V_c_init; T_init];

% Solve the system using ode45
[t, y] = ode45(ode_system, t_span2, initial_conditions);

% Extract the solutions for V_c and T
V_c_2 = y(:, 1); % Voltage across the capacitor
T_2 = y(:, 2);   % Temperature



figure;
subplot(1,2,1)
plot(t, V_c_2, 'b-', 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Voltage across capacitor V_c(t)');
grid on;

subplot(1,2,2)
plot(t, T_2, 'r-', 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Temperature T(t)');
title('Temperature vs. Time');
grid on;


%% DC I-V characteristics
%define parameters
g_0 = 5e-3;
g_1 = 1700;
C_th = 1e-14;
g_th= 1/1.5e6;
T_0 = 300;
V_b = 1.1;
R_b = 1e3;
C = 10e-9;
R_s = 200;


% Temperature swept
T = [300.01 : 0.01 : 7000];

R_s = 0;

% Set dT/dt = 0
vmc = sqrt(g_th .* (T-T_0) ./ G(T));
im_R0 = vmc .* G(T);
vm_R0 = vmc .* (1 + R_s .* G(T))

R_s = 200;

im_R200 = vmc .* G(T);
vm_R200 = vmc .* (1 + R_s .* G(T))

figure;
plot(vm_R0, im_R0, 'b-', 'LineWidth', 2);
hold on;
plot(vm_R200, im_R200, 'r-', 'LineWidth', 2);
hold off;
xlabel('Im(ma)');
ylabel('Vm(V)');
grid on;


%%Single Cell with Rs
%define parameters
g_0 = 5e-3;
g_1 = 1700;
C_th = 1e-14;
g_th= 1/1.5e6;
T_0 = 300;
V_b = 1.1;
R_b = 1e3;
C = 10e-9;
R_s = 200;
R_c = 10;
%Initial conditions
T_init = 300;
V_c_init = 0.001;

%Timespan
t_span3 = [0 : 0.00001e-4 : 2.5e-4];

%Functions
G = @(T) g_0 .* exp(-g_1 ./ T);
V_mc = @(V_c, T) V_c ./ (1 + G(T) .* R_s);
%Define system of ODEs
ode_system = @(t, y) [(((V_b - y(1)) ./ R_b - G(y(2)) .* y(1)) - (y(1)-y(3))./R_c)./ C; 
                      ((V_mc(y(1), y(2)) .^ 2) .* G(y(2)) - (y(2)-T_0) .* g_th) ./ C_th;
                      (((V_b - y(3)) ./ R_b - G(y(4)) .* y(3)) - (y(3)-y(1))./R_c)./ C; 
                      ((V_mc(y(3), y(4)) .^ 2) .* G(y(4)) - (y(4)-T_0) .* g_th) ./ C_th
                      ];
% Initial conditions vector [V_c, T]
initial_conditions = [V_c_init; T_init; V_c_init; T_init];

% Solve the system using ode45
[t, y] = ode45(ode_system, t_span3, initial_conditions);

% Extract the solutions for V_c and T
V_c_3 = y(:, 1); % Voltage across the capacitor
T_3 = y(:, 2);   % Temperature
V_c_4 = y(:, 3); % Voltage across the capacitor
T_4 = y(:, 4);   % Temperature



figure;
subplot(1,2,1)
plot(t, V_c_3, 'b-', 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Voltage across capacitor V_c(t)');
grid on;

subplot(1,2,2)
plot(t, V_c_4, 'r-', 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Voltage across capacitor V_c(t)');
grid on;