
close all
clear


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
R_c = 10000;
%Initial conditions
T_init = 300;
n_iterations = 5;
V_c_init2 = 0.68; %sweep init voltage of capacitor in cell 2 from 0 to 1V in 5 steps
colorstring = 'rgbcm';
V_c_init1_const = 0.0; %init voltage of capacitor cell 1
sim_time = 70e-6;
%Timespan
t_span = [0 : 0.000007e-6 : sim_time];

%Functions
G = @(T) g_0 .* exp(-g_1 ./ T);
V_mc = @(V_c, T) V_c ./ (1 + G(T) .* R_s);


figure;
hold on;
title('Sweep over Vc init');
%subplot(1,2,1);
xlabel('Time in seconds');
ylabel('Voltage across capacitor V_c(t)');
grid on;


    %Define system of ODEs
    ode_system = @(t, y) [(((V_b - y(1)) ./ R_b - G(y(2)) .* y(1)) - (y(1)-y(3))./R_c)./ C; 
                          ((V_mc(y(1), y(2)) .^ 2) .* G(y(2)) - (y(2)-T_0) .* g_th) ./ C_th;
                          (((V_b - y(3)) ./ R_b - G(y(4)) .* y(3)) - (y(3)-y(1))./R_c)./ C; 
                          ((V_mc(y(3), y(4)) .^ 2) .* G(y(4)) - (y(4)-T_0) .* g_th) ./ C_th
                          ];
    % Initial conditions vector [V_c, T]
    initial_conditions = [V_c_init1_const; T_init; V_c_init2; T_init]; %sweep V_c_init for cell 1, V_c_init for cell 2 remains constant
    
    % Solve the system using ode45
    [t, y] = ode45(ode_system, t_span, initial_conditions);
    
    % Extract the solutions for V_c and T
    V_c_1 = y(:, 1); % Voltage across the capacitor in cell 1 
    T_1 = y(:, 2);   % Temperature
    V_c_2 = y(:, 3); % Voltage across the capacitor in cell 2
    T_2 = y(:, 4);   % Temperature

    % Plot the result
    plot(t, V_c_1, 'r-', 'LineWidth', 1)
    plot(t, V_c_2, 'm-', 'LineWidth', 1)


    %%Calculate Minima for V_c_1
    T = 1.416e-5 %read from graph! Not calculated


    %V_c_1 has 10000001 element
    idx1 = round((4/7)*length(V_c_1));
    %idx2 = round(5/7*length(V_c_1));
    idx2 = idx1 + round((T/sim_time) * length(V_c_1));
    [sectionMin, sectionRelIdx] = min(V_c_1(idx1:idx2));  % Minimum 
    %sectionAbsIdx = sectionRelIdx + 4/7*length(V_c_1) - 1;  % Convert to absolute index
    %disp(sectionMin)
    absIdx_V1 = idx1 + sectionRelIdx;
    abstime_V1 = (absIdx_V1/length(V_c_1)) * sim_time;
    absminvalue_V1 = sectionMin;

    %%Calculate Minima for V_c_2

    idx1 = round(4/7*length(V_c_2));
    %idx2 = round(5/7*length(V_c_2));
    idx2 = idx1 + round((T/sim_time) * length(V_c_2));
    [sectionMin, sectionRelIdx] = min(V_c_2(idx1:idx2));  % Minimum 
    %sectionAbsIdx = sectionRelIdx + 4/7*length(V_c_1) - 1;  % Convert to absolute index
    %disp(sectionMin)
    absIdx_V2 = idx1 + sectionRelIdx;
    abstime_V2 = (absIdx_V2/length(V_c_2)) * sim_time;
    absminvalue_V2 = sectionMin;
    
    

    phaseshift = ((abstime_V2 - abstime_V1)/T) *360
