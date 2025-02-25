tic;

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

%% V_c and resistor values for different input images

%V_c
v_data = load('V_c_init_col_mat.mat');
input_v = v_data.V_c_init_col_mat;
image_no = 2; %from 1 to 7

var_mat = extract_vc(input_v, image_no);




%Resistor values

%Mapping function: Exponential
% r_data  = load('exp_weights_mat.mat');
% inp = r_data.exp_res_mat;
% 
%Mapping function: Sigmoid
% r_data  = load('sigmoid_weights_mat.mat');
% inp = r_data.sig_res_mat;

% Mapping function: Linear
% r_data  = load('linear_res_mat.mat');
% inp = r_data.linear_res_mat;
% 
% Mapping function: Tuned Linear
% r_data  = load('tune_linear_res_mat.mat');
% inp = r_data.tune_linear_res_mat;

% %Mapping function: Medium Linear
% r_data  = load('med_linear_res_mat.mat');
% inp = r_data.med_linear_res_mat;

%Mapping function: Medium LInear Test
r_data = load("med_linear_res_mat_test.mat")
inp = r_data.med_linear_res_mat_test;

n = size(inp, 1);
inp(1:n+1:end) = 1;

res_mat = inp;
%imagesc(res_mat)
%% Random values for V_C and resistor matrix
% 
% %generating random resistance and vc values for n cells
%n = 49; % no of cells
% % 
% % A = rand(n);
% % A = 10 + (A * (50 - 10));
% % A(1:10+1:end) = 1;
% % A = (A + A') / 2;
% % res_mat = A .* 1000;
% % 
% B = [rand(1, n); 300*ones(1,n)];
% 
% 
% a = 0.4; % First alternate value
% b = 0.7; % Second alternate value
% 
% B(1, 1:2:end) = a * rand(1, ceil(size(B, 2) / 2)); % Set elements at odd indices to 'a'
% B(1, 2:2:end) = b + (1 - b) * rand(1, floor(size(B, 2) / 2)); % Set elements at even indices to 'b'
% 
% var_mat = B;

%imagesc(res_mat)
%imagesc(var_mat(1,:))
%% User-defined values for V_c and resistor matrix

% %these are the values for the resistor between two cells
% res_mat =       [ 1, 2.5;
%                   2.5, 1] .* 1000; 
% % 
% % 
% 
% 
% 
% %Initial conditions
% var_mat    =    [ 0, 300;     % V_c1, T_1
%                  0.69, 300]';    

%% Running the ODE for the network

var_flat = var_mat(:); %flatting the intial condition matrix for input in ODE
t_span = [0:0.0001e-6:2.5e-4]; %time span of simulation

odeFunc = @(t, y) odeMatrix(t, y, res_mat, params);
%[t, y] = ode45(odeFunc, t_span, var_flat); %ode45 solver
[t, y] = ode15s(odeFunc, t_span, var_flat); %ode15s solver
if isempty(y)
    error('Simulation returned no data.');
end

Y = reshape(y', 2, [], length(t)); %reshape output to a 2*N matrix

disp("ODEs solved:")
toc;
disp(' ')

%% Plotting graphs

%plot_sep(Y, t); %plot separate graphs
plot_overlap(Y, t); %plot overlapping graphs

disp("Plotting completed:")
toc;
disp(' ')

%% Calculating attributes (such as amplitude, time period, frequency, phase shifts, pixel values)

attr_mat = attrib(Y, t);
%disp(attr_mat);



disp("Calculation of attributes completed:")
toc;
disp(' ')

%% Displaying the attributes as a look-up table

%matrixToTable(attr_mat, var_mat)

%% Plotting pixel values as images

% %plotting a single image
% pixel_mat = attr_mat(5,:);
% image(pixel_mat)

%plotting images over time (multiple)
multi_image(attr_mat)

disp("Plotting multiple images completed:")
toc;
disp(' ')

%% FUNCTION definitions

%Function to convert attributes matrix to a table
function matrixToTable(attributes, var_mat)

    % Extract the first row from matrix a
    rowFromA = var_mat(1, :);
    
    % Extract the 1st, 2nd, 3rd, 22nd, and 23rd rows from matrix b
    rowsFromB = attributes([1, 2, 3, 22, 23], :);
    
    % Combine the rows (row from a followed by rows from b)
    combinedMatrix = [rowFromA; rowsFromB];

    rowNames = {'Inital V_c', 'Amplitudes', 'Time Period', 'Frequency', 'Phase shift wrt Cell 1', 'Pixel Values'};
    
    % Define column names dynamically as 'cell 1', 'cell 2', ..., 'cell N'
    numCols = size(combinedMatrix, 2);
    colNames = strcat('Cell ', string(1:numCols));
    
    % Convert the matrix to a table with specified row and column names
    table = array2table(combinedMatrix, 'RowNames', rowNames, 'VariableNames', colNames);
    % Display the table
    disp(table);
end

%Function to extract V_C values for a particular image (used in testing)
function vc_mat = extract_vc(inputMatrix, image_no)
    vc_mat(1, :) = [ inputMatrix(:,image_no) ]';
    vc_mat(2, :) = 300*ones(1,size(inputMatrix, 1));
end

%function to plot a multiple image
function multi_image(inputMatrix)
    % This function extracts rows 5, 7, 9, 11, ..., 23 from the input matrix
    rowsToExtract = 5:2:23; % Generate the row indices (5, 7, 9, ..., 23)
    
    % Extract the corresponding rows from the input matrix
    pixel_mat = inputMatrix(rowsToExtract, :);
    
    
    % Line thickness (in pixels) - set this to a very small value
    lineThickness = 1; % Thin line thickness
    
    % Set the invisible line color to a custom RGB color (light gray here)
    invisibleLineColor = [0.9, 0.9, 0.9]; % Light gray RGB value (you can set any color)
    
    % Initialize the final concatenated image
    finalImage = [];
    
    % Loop through each row of the matrix
    for i = 1:size(pixel_mat, 1)
        % Step 1: Reshape the current row into a 7x7 matrix
        image_mat = reshape(pixel_mat(i, :), 7,7);
        
        % Step 2: Convert the binary image to RGB
        rgbImage = cat(3, image_mat, image_mat, image_mat); % Convert to grayscale RGB
        
        % Step 3: Add the current image to the final image with an invisible line
        if isempty(finalImage)
            finalImage = rgbImage; % Initialize with the first image
        else
            % Create a custom colored line (invisible in light gray)
            invisibleLine = repmat(reshape(invisibleLineColor, 1, 1, 3), size(rgbImage, 1), lineThickness); 
            
            % Concatenate the invisible line and the next image
            finalImage = [finalImage, invisibleLine, rgbImage]; % Concatenate horizontally
        end
    end
    
    % Display the final concatenated image with custom invisible lines
    figure;
    imshow(finalImage, 'InitialMagnification', 'fit'); % Scaled to fit screen
    %set(gcf, 'Color', [1 1 1]); % Set figure background color to white (default)
    title('Image Correction over Time');
    % 
    % % Save the final concatenated image as a PNG
    % imwrite(finalImage, 'concatenatedImageWithCustomInvisibleLines.png');
    % 
    % disp('Final concatenated image with custom invisible lines saved as concatenatedImageWithCustomInvisibleLines.png');
end

%function to plot a single image
function image(pixel_mat)
    figure();
    %convert pixel values to images
    image_mat = reshape(pixel_mat, 6,5)
    
    % Step 2: Convert to logical format (optional, for clarity)
    binaryImage = logical(image_mat);
    
    % Step 3: Display the binary image
    imshow(binaryImage, 'InitialMagnification', 5000); % Scaled up for visibility
    title('Binary Image: Black and White');
    
    % Step 4: Save the image as a binary PNG
    imwrite(binaryImage, 'binaryImage.png');
    
    disp('Image saved as binaryImage.png');
end

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
    
    attr_mat = zeros(23, size(inp_mat, 1));
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
        attr_mat(2, a) = mean(diff(locs)) * 1e-10; % Average time difference between peaks
    
        % frequency
        attr_mat(3, a) = 1 ./ attr_mat(2, a);
        %disp(attr_mat);
          
    end
    
    %for phase diff
    t1= 1.5e-5;
    av_tp = sum(attr_mat(2, :)) / size(attr_mat, 2);
    
    
    for num = 1:10
        
        t2 = t1 + av_tp;
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

%             disp("current loop")
%             disp(a)
            

            if length(loc1) == length(loc2)
                delta = loc2 - loc1;
            else
                disp('No of peaks of wave 1');
                disp(num2str(length(loc1)));
                disp('No of peaks of wave 2');
                disp(num2str(length(loc2)));
                
            end
    
            phase_diff = 360 * delta / av_tp;
            
            phase_scaled = mod(phase_diff + 180, 360) - 180; %phase difference in [-180, 180]
            %phase_scaled = mod(phase_diff, 360); %phase difference in [0, 360]
            
            r1 = 2 * num +2;
            r2= 2* num + 3;
            attr_mat(r1, a) = phase_scaled;
    
            %for pixel values
            if phase_scaled <= 0 
                attr_matr(r2, a) = 0;    %black pixel
            else
                attr_mat(r2, a) = 1;     %white pixel
            end
    
        end
        %disp(attr_mat);
        t1 = t2;
        
    end

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
        %ylim([0.8 1]);
        grid on;
        hold on;

    end
end