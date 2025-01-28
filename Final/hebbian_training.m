%% 49x49 Weights Training Matrix Calculation 

clear

% Define image dimensions and number of training images
rows = 7;
cols = 7;
num_images = 5;
num_pixels = rows * cols;

% Initialize a matrix to store the flattened training images
training_images = zeros(num_pixels, num_images);

filename = {'number0.txt', 'number1.txt', 'number3.txt', 'number4.txt', 'number7.txt'};


% Load 5 binary images
for i = 1:5
    data = readmatrix(filename{i});
    disp('data is:');
    disp(data);
    img = imbinarize(data);  % Ensure binary format
    img_vector = img(:);    % Flatten the image into a vector
    training_images(:, i) = img_vector;
end

% Initialize the weights matrix
weights_matrix = zeros(num_pixels, num_pixels);

% Compute the weights matrix using the Hebbian Learning Rule
for i = 1:num_pixels
    for j = 1:num_pixels
        if i ~= j
            % Compute the weight between pixel i and pixel j
            for col=1:5
                weights_matrix(i, j) = weights_matrix(i, j) + training_images(i, col) .* training_images(j, col);
            end
        end
    end
end

weights_matrix = weights_matrix./num_images;

% Display the resulting weights matrix
disp('Weights matrix:');
disp(weights_matrix);

% Save the weights matrix to a file
save('weights_matrix.mat', 'weights_matrix');

imagesc(weights_matrix);


%% Extracting the Weights Vector (Upper Triangle of the Weights Matrix) %%

% Extract the upper right triangle (excluding the diagonal)
upper_triangle_indices = triu(true(num_pixels), 1); % Create a logical mask for upper triangle
weights_vector = weights_matrix(upper_triangle_indices); % Extract the values

% Display the vector of weights
disp('Upper right triangle weights vector:');
disp(weights_vector);

% Save the weights vector to a file
save('weights_vector.mat', 'weights_vector');

wmax = max(weights_vector);
wmin = min(weights_vector);
fprintf('The Min and Max of weights are: (%d, %d)\n', wmin, wmax);


%% Map to Resistors %%

% Linear Mapping %%%%%%%%%%%%%%%%%%%%%
linear_res = weights_vector.*100000+10000;
linear_res_mat = weights_matrix.*100000+10000;
minValue = min(linear_res(:)); % Overall minimum
maxValue = max(linear_res(:)); % Overall maximum

fprintf('The Min and Max of normal Linear are: (%d, %d)\n', minValue, maxValue);

% Save the Linear-transformed weights matrix to a file
save('linear_res_mat.mat', 'linear_res_mat');

% Plot the Color Map
figure(1);
imagesc(linear_res_mat);
hold on;

% Plot the Linear Transform Graph
figure(2);
plot(weights_vector,linear_res);
hold off;


%% Sigmoid Mappin %%%%%%%%%%%%%%%%%%%%
sig_res =  (1./(1+exp(-(-6+12.*weights_vector)))).*100000;
sig_res_mat = (1./(1+exp(-(-6+12.*weights_matrix)))).*100000;
minValue2 = min(sig_res(:)); % Overall minimum
maxValue2 = max(sig_res(:)); % Overall maximum

fprintf('The Min and Max of Sigmoid are: (%d, %d)\n', minValue2, maxValue2);

% Save the sigmoid-transformed weights matrix to a file
save('sigmoid_weights_mat.mat', 'sig_res_mat');

% Plot the Color Map
figure(1);
imagesc(sig_res_mat);
hold on;

% Plot the Linear Transform Graph
figure(2);
plot(weights_vector,sig_res);
hold off;


%% Exponential Mappin %%%%%%%%%%%%%%%%%%%%
exp_res =  (exp(weights_vector.*10)).*100;
exp_res_mat = (exp(weights_matrix.*10)).*100;
minValue3 = min(exp_res(:)); % Overall minimum
maxValue3 = max(exp_res(:)); % Overall maximum

fprintf('The Min and Max of Exponential are: (%d, %d)\n', minValue3, maxValue3);

% Save the exp-transformed weights matrix to a file
save('exp_weights_mat.mat', 'exp_res_mat');

% Plot the Color Map
figure(1);
imagesc(exp_res_mat);
hold on;

% Plot the exp Transform Graph
figure(2);
plot(weights_vector,exp_res);
hold off;


%% Parabular Mapping %%%%%%%%%%%%%%%%%%%%%
par_res = (weights_vector.^4).*100000 + weights_vector.*10000 + 4560;
par_res_mat = (weights_matrix.^4).*100000 + weights_matrix.*10000 + 4560;
minValue6 = min(par_res(:)); % Overall minimum
maxValue6 = max(par_res(:)); % Overall maximum

fprintf('The Min and Max of Parabular are: (%d, %d)\n', minValue6, maxValue6);

% Save the parabular-transformed weights matrix to a file
save('par_res_mat.mat', 'par_res_mat');

% Plot the Color Map
figure(1);
imagesc(par_res_mat);
hold on;

% Plot the exp Transform Graph
figure(2);
plot(weights_vector,par_res);
hold off;


%% Logarithmic Mapping %%%%%%%%%%%%%%%%%%%%%
log_res = log(10+weights_vector.*100).*12345;
log_mat = log(10+weights_matrix.*100).*12345;
minValue7 = min(log_res(:)); % Overall minimum
maxValue7 = max(log_res(:)); % Overall maximum

fprintf('The Min and Max of Logarithmic are: (%d, %d)\n', minValue7, maxValue7);

% Save the Log-transformed weights Matrix to a file
save('log_mat.mat', 'log_mat');

% Plot the Color Map
figure(1);
imagesc(log_mat);
hold on;

% Plot the exp Transform Graph
figure(2);
plot(weights_vector,log_res);
hold off;

