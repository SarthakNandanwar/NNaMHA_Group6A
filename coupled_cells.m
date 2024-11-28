clc;
clear;
close all;



%define parameters under a struct (so it can be accessed in any functions by passing the struct as an input)
params = struct();

params.g_0  = 5e-3;
params.g_1  = 1700;
params.C_th = 1e-14;
params.g_th = 1/1.5e6;
params.T_0  = 300;
params.V_b  = 1.1;
params.R_b  = 1e3;
params.C    = 10e-9;
params.R_s  = 200;

params.res_mat =   [1, 2, 3, 2, 5;
                    2, 1, 3, 4, 5;
                    3, 3, 1, 3, 5;
                    2, 4, 3, 1, 5;
                    5, 5, 5, 5, 1] .* 100;


%Initial conditions
var_mat    =    [0.2, 300;      % V_c1, T_1
                 0.3, 300;      % V_c2, T_2
                 0.1, 300;      % V_c3, T_3
                 0.4, 300;      % V_c4, T_4
                 0.5, 300]';    % V_c5, T_5

var_flat = var_mat(:); %flatting the intial condition matrix for input in ODE

t_span = [0:0.0001e-4:2.5e-4];


odeFunc = @(t, y) odeMatrix(t, y, params);
[t, y] = ode45(odeFunc, t_span, var_flat);

Y = reshape(y', 2, [], length(t)); %reshape output to a 2*N matrix


%plotting separate graphs for each circuit
figure;
title('Solutions of the Differential Equations');

for i = 1:size(var_mat,2)
    subplot(size(var_mat,2), 1, i);
    plot(t, squeeze(Y(1,i,:)));
    
    xlabel('Time');
    ylabel('V_'+string(i));
    title('Values of V_'+string(i)');
    ylim([0.8 1]);
    grid on;
    %hold on;
end

%plotting overlapping graphs for each circuit
figure;
title('Solutions of the Differential Equations');

for i = 1:size(var_mat,2)
    plot(t, squeeze(Y(1,i,:)));
    
    xlabel('Time');
    ylabel('V');
    title('Values of V');
    %xlim([0 0.2]);
    %ylim([0.8 1]);
    grid on;
    hold on;
end

%ode function
function dydt = odeMatrix(t, y, params)

    %define a new matrix that converts y into a matrix
    Y_mat = reshape(y, 2, []);
    output = zeros(size(Y_mat));


    sub = calc_sub(Y_mat, params);
    G = @(T, params) params.g_0 * exp(-params.g_1 / T);
    V_mc = @(V_c, T, params) V_c / (1 + G(T, params) * params.R_s);


    %connect pairs
    for i = 1:size(Y_mat,2)
        V_c = Y_mat(1,i);
        T = Y_mat(2,i);

        output(:, i) = [  ( ((params.V_b - V_c)/params.R_b)  - (G(T, params) * V_mc(V_c, T, params)) - sub(1,i) ) / params.C;
                          ( ((V_mc(V_c, T, params) ^ 2) * G(T, params))  - ((T - params.T_0) * params.g_th )) / params.C_th];
    end
    
    dydt = output(:);
end

function out_mat = calc_sub(var_mat, params)
    
    out_mat = zeros(1, size(var_mat,2));

    for i = 1:size(var_mat, 2)
        for j = 1:size(var_mat, 1)
            out_mat(1, i) = out_mat(1, i) + ( (var_mat(1, i) - var_mat(1, j)) / params.res_mat(i,j) ) ;
        end
    end

end