clc;
clear;
close all;

format long g;

syms x;

digits(20);

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
%params.t_delay = 0.1e-4;

%generating random resistance and vc values for n cells
n = 20; % no of cells

A = rand(n);
A = 60 + (A * (80 - 60));
A(1:10+1:end) = 1;
A = (A + A') / 2;
res_mat = A .* 1000;

B = [rand(1, n); 300*ones(1,n)];
var_mat = B;

% %these are the values for the resistor between two cells
% res_mat =       [ 1, 80, 50, 65;      %resistance values for connections to the first cell
%                  80,  1, 70, 78;      %resistance values for connections to the second cell and so on
%                  50, 70,  1, 30;
%                  65, 78, 30,  1;] .* 1000; 
% 
% 
% 
% 
% 
% %Initial conditions
% var_mat    =    [ 0.1, 300;     % V_c1, T_1
%                  0.44, 300;     % V_c2, T_2 and so on
%                  0.67, 300;
%                  0.98, 300]';    



var_flat = var_mat(:); %flatting the intial condition matrix for input in ODE

%initial time delays for inputs to each circuit
%t_delay = [0, 0, 0];

t_span = [0:0.0001e-4:2.5e-4]; %time span of simulation


odeFunc = @(t, y) odeMatrix(t, y, res_mat, params);
[t, y] = ode45(odeFunc, t_span, var_flat);

Y = reshape(y', 2, [], length(t)); %reshape output to a 2*N matrix
%disp(size(Y, 2));

attr_mat = attrib(Y, t);
%disp(attr_mat);

%for displaying the matrix with labels
rows = {'Amplitude(V)', 'Time Period(s)', 'Frequency(Hz)', 'Phase Shift wrt V1(deg)', 'Pixel Value'};
cols = arrayfun(@(x) ['CELL ', num2str(x)], 1:size(attr_mat, 2), 'UniformOutput', false);
attr_table = array2table(attr_mat, 'VariableNames', cols, 'RowNames', rows);
disp(attr_table);


%plot_sep(Y, t); %plot separate graphs
plot_overlap(Y, t); %plot overlapping graphs


%ode function
function dydt = odeMatrix(t, y, res_mat, params)

    %define a new matrix that converts y into a matrix
    Y_mat = reshape(y, 2, []);
    output = zeros(size(Y_mat));
    %V_b = zeros(size(time_delay));


    sub = calc_sub(Y_mat, res_mat);
    G = @(T, params) params.g_0 * exp(-params.g_1 / T);
    V_mc = @(V_c, T, params) V_c / (1 + G(T, params) * params.R_s);
    %V_b = @(t) 1.1 * heaviside(t - 0.1e-4);


    %connect pairs
    for i = 1:size(Y_mat,2)
        V_c = Y_mat(1,i);
        T = Y_mat(2,i);
        %V_b = @(t) (t >= t_delay(i)) * 1.1;
        %disp(V_b);

        output(:, i) = [  ( ((params.V_b - V_c)/params.R_b)  - (G(T, params) * V_mc(V_c, T, params)) - sub(1,i) ) / params.C;
                          ( ((V_mc(V_c, T, params) ^ 2) * G(T, params))  - ((T - params.T_0) * params.g_th )) / params.C_th];

%           output(:, i) = [  ( ((V_b(t) - V_c)/params.R_b)  - (G(T, params) * V_mc(V_c, T, params)) - sub(1,i) ) / params.C;
%                            ( ((V_mc(V_c, T, params) ^ 2) * G(T, params))  - ((T - params.T_0) * params.g_th )) / params.C_th];
    end
    
    dydt = output(:);
end

function out_mat = calc_sub(var_mat, res_mat)
    
    out_mat = zeros(1, size(var_mat,2));

    for i = 1:size(var_mat, 2)
        for j = 1:size(var_mat, 1)
            out_mat(1, i) = out_mat(1, i) + ( (var_mat(1, i) - var_mat(1, j)) / res_mat(i,j) ) ;
        end
    end
    

end

%calculating the attributes(amplitude, time period, frequency, and PHASE SHIFT) for different signals
function attr_mat = attrib(inp_mat, t)
    
    attr_mat = zeros(5, size(inp_mat, 1));
    n = size(inp_mat, 2); % no of cells

    t_0 = 0.0;  % Start looking for the lag after 2 seconds
    start_idx = find(t >= t_0, 1);  % Find the index where time is greater than or equal to t_0
    
    for a = 1:n

        y = squeeze(inp_mat(1,a,:));
        %y = inp_mat(a,:);
    
        % Find the peak amplitude
        attr_mat(1, a) = max(y);
    
        % Find the time period
        [pks, locs] = findpeaks(y, start_idx); % Find peaks and their time indices
        attr_mat(2, a) = mean(diff(locs)) * 1e-8; % Average time difference between peaks
    
        % frequency
        attr_mat(3, a) = 1 ./ attr_mat(2, a);
        %disp(attr_mat);
          
    end
    
    %for phase diff
    
    t1= 2.14e-4;
    av_tp = sum(attr_mat(2, :)) / size(attr_mat, 2);
    t2 = t1+ av_tp;

    t_ind = (t>= t1) & (t <=t2);
    t_rest = t(t_ind);

    for a = 1:n
%         ref_wave = inp_mat(1, :);
%         curr_wave = inp_mat(a, :);
        
        ref_wave = squeeze(inp_mat(1,1,:));
        curr_wave =squeeze(inp_mat(1,a,:));

        ref_rest = ref_wave(t_ind);
        curr_rest = curr_wave(t_ind);

        [peak1, loc1] = findpeaks(ref_rest, t_rest);
        [peak2, loc2] = findpeaks(curr_rest, t_rest);
        
        if length(loc1) == length(loc2)
            delta = loc2 - loc1;
        else
            disp('No of peaks of wave 1');
            disp(num2str(length(loc1)));
            disp('No of peaks of wave 2');
            disp(num2str(length(loc2)));
            
        end

        phase_diff = 360 * delta / av_tp;
        
        %phase_scaled = mod(phase_diff + 180, 360) - 180; %phase difference in [-180, 180]
        phase_scaled = mod(phase_diff, 360); %phase difference in [0, 360]

        attr_mat(4, a) = phase_scaled;

        %for pixel values
        if phase_scaled >= 180 
            attr_matr(5, a) = 0;    %black pixel
        else
            attr_mat(5, a) = 1;     %white pixel
        end

    end
    %disp(attr_mat);
   

end

%function for plotting separate graphs for each circuit
function plot_sep(Y, t)
    figure('Position', [20, 270, 700, 500]);
    title('Solutions of the Differential Equations');
    
    for i = 1:size(Y,2)
        subplot(size(Y,2), 1, i);
    
        curr_vals = squeeze(Y(1,i,:));
        
        plot(t, curr_vals);
        
    %     [max_val, max_time] = findpeaks(squeeze(Y(1,i,:)), t)
    %     plot(max_time, max_val, '^r', 'MarkerFaceColor','g');
    
        xlabel('Time');
        ylabel('V_'+string(i));
        title('Values of V_'+string(i)');
        %ylim([0.8 1]);
        grid on;
        %hold on;
    end
end

%plotting overlapping graphs for each circuit
function plot_overlap(Y, t)
    
    figure('Position', [50, 270, 1500, 300]);
    title('Solutions of the Differential Equations');
    
    for i = 1:size(Y,2)

        curr_vals = squeeze(Y(1,i,:));
        plot(t, curr_vals);
    %     [max_val, max_time] = findpeaks(squeeze(Y(1,i,:)), t)
    %     plot(max_time, max_val, '^r', 'MarkerFaceColor','g');
        xlabel('Time');
        ylabel('V');
        title('Values of V');
        %xlim([0 0.2]);
        ylim([0.8 1]);
        grid on;
        hold on;

    end
end