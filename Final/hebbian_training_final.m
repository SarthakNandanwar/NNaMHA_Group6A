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
    img = imbinarize(data);  % Ensure binary format (optional)
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

weights_matrix = weights_matrix/num_pixels;

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

% Save the weights matrix and vector to a file
% save('weights_matrix.mat', 'weights_matrix');
save('weights_vector.mat', 'weights_vector'); % Save the vector

%% Map to Resistors %%

% Linear Mapping %%%%%%%%%%%%%%%%%%%%%
linear_res = weights_vector.*100000+10000;
linear_res_mat = weights_matrix.*100000+10000;
minValue = min(linear_res(:)); % Overall minimum
maxValue = max(linear_res(:)); % Overall maximum

fprintf('The Min and Max of normal Linear are: (%d, %d)\n', minValue, maxValue);

% disp(linear_res)
%fprintf('Resistors Linear Mapping of Weights are: \n [%s]\n', join(string(linear_res_mat), ', '));

% Save the Linear-transformed weights matrix to a file
save('linear_res_mat.mat', 'linear_res_mat');
imagesc(linear_res_mat);



%% Sigmoid Mappin %%%%%%%%%%%%%%%%%%%%
sig_res =  (1 ./ (1 + exp(-weights_vector)))*10000;
sig_res_mat =  (1 ./ (1 + exp(-weights_matrix)))*10000;
minValue2 = min(sig_res(:)); % Overall minimum
maxValue2 = max(sig_res(:)); % Overall maximum

fprintf('The Min and Max of Sigmoid are: (%d, %d)\n', minValue2, maxValue2);

% Display the resulting sigmoid-transformed weights matrix
%fprintf('Resistors Sigmoid Mapping of Weights are: \n [%s]\n', join(string(sig_res_mat), ', '));

% Save the sigmoid-transformed weights matrix to a file
save('sigmoid_weights_mat.mat', 'sig_res_mat');
imagesc(sig_res_mat);


%% Exponential Mappin %%%%%%%%%%%%%%%%%%%%
exp_res =  10000*exp(weights_vector.*10);
exp_res_mat =  10000*exp(weights_matrix.*10);
minValue3 = min(exp_res(:)); % Overall minimum
maxValue3 = max(exp_res(:)); % Overall maximum

fprintf('The Min and Max of Exponential are: (%d, %d)\n', minValue3, maxValue3);

% Display the resulting sigmoid-transformed weights matrix
%fprintf('Resistors Exponential Mapping of Weights are: \n [%s]\n', join(string(exp_res_mat), ', '));

% Save the sigmoid-transformed weights matrix to a file
save('exp_weights_mat.mat', 'exp_res_mat');
imagesc(exp_res_mat);



%% Tuned Linear Mapping %%%%%%%%%%%%%%%%%%%%%
tune_linear_res = (weights_vector.*1000000 + 10000)*0.8;
tune_linear_res_mat = (weights_matrix.*1000000 + 10000)*0.8;
minValue4 = min(tune_linear_res(:)); % Overall minimum
maxValue4 = max(tune_linear_res(:)); % Overall maximum

fprintf('The Min and Max of Tuned Linear are: (%d, %d)\n', minValue4, maxValue4);

%disp(linear_res)
%fprintf('Resistors Tuned Linear Mapping of Weights are: \n [%s]\n', join(string(tune_linear_res_mat), ', '));

% Save the Linear-transformed weights matrix to a file
% save('tune_linear_res.mat', 'tune_linear_res');
save('tune_linear_res_mat.mat', 'tune_linear_res_mat');
imagesc(tune_linear_res_mat);


%% Medium Linear Mapping %%%%%%%%%%%%%%%%%%%%%
med_linear_res = (weights_vector.*100000 + 10000)*2.4;
med_linear_res_mat = (weights_matrix.*100000 + 10000)*2.4;
minValue5 = min(med_linear_res(:)); % Overall minimum
maxValue5 = max(med_linear_res(:)); % Overall maximum

fprintf('The Min and Max of Medium Linear are: (%d, %d)\n', minValue5, maxValue5);

% disp(linear_res)
% fprintf('Resistors Medium Linear Mapping of Weights are: \n [%s]\n', join(string(med_linear_res), ', '));
%fprintf('Resistors Matrix Medium Linear Mapping of Weights are: \n [%s]\n', join(string(med_linear_res_mat), ', '));

% Optional: Save the Linear-transformed weights vector to a file
% save('med_linear_res.mat', 'med_linear_res');
save('med_linear_res_mat.mat', 'med_linear_res_mat');
imagesc(med_linear_res_mat);

%% Parabular Mapping %%%%%%%%%%%%%%%%%%%%%
par_res = (weights_vector.^2).*100000 + weights_vector.*10000 + 43560;
par_res_mat = (weights_matrix.^2).*100000 + weights_matrix.*10000 + 43560;
minValue6 = min(par_res(:)); % Overall minimum
maxValue6 = max(par_res(:)); % Overall maximum

fprintf('The Min and Max of Parabular are: (%d, %d)\n', minValue6, maxValue6);

% disp(linear_res)
% fprintf('Resistors Medium Linear Mapping of Weights are: \n [%s]\n', join(string(med_linear_res), ', '));
%fprintf('Resistors Matrix Medium Linear Mapping of Weights are: \n [%s]\n', join(string(med_linear_res_mat), ', '));

% Optional: Save the Linear-transformed weights vector to a file
% save('med_linear_res.mat', 'med_linear_res');
save('par_res_mat.mat', 'par_res_mat');
imagesc(par_res_mat);

%% Logarithmic Mapping %%%%%%%%%%%%%%%%%%%%%
log_res = 12345*log((weights_vector.^2).*100000 + weights_vector.*10000 + 43560 );
log_mat = 12345*log((weights_matrix.^2).*100000 + weights_matrix.*10000 + 43560);
minValue7 = min(log_res(:)); % Overall minimum
maxValue7 = max(log_res(:)); % Overall maximum

fprintf('The Min and Max of Logarithmic are: (%d, %d)\n', minValue7, maxValue7);

% disp(linear_res)
% fprintf('Resistors Medium Linear Mapping of Weights are: \n [%s]\n', join(string(med_linear_res), ', '));
%fprintf('Resistors Matrix Medium Linear Mapping of Weights are: \n [%s]\n', join(string(med_linear_res_mat), ', '));

% Optional: Save the Linear-transformed weights vector to a file
% save('med_linear_res.mat', 'med_linear_res');
save('log_mat.mat', 'log_mat');
imagesc(log_mat);


% %% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clear
% 
% % Define image dimensions and number of training images
% rows = 2;
% cols = 2;
% num_images = 3;
% num_pixels = rows * cols;
% 
% % Initialize a matrix to store the flattened training images
% training_images = zeros(num_pixels, num_images);
% 
% filename = {'Test1.txt', 'Test2.txt', 'Test3.txt'};
% 
% 
% % Load 3 binary images
% for i = 1:3
%     data = readmatrix(filename{i});
%     disp('data is:');
%     disp(data);
%     img = imbinarize(data);  % Ensure binary format (optional)
%     img_vector = img(:);    % Flatten the image into a vector
%     disp(img_vector);
%     training_images(:, i) = img_vector;
% end
% 
% disp('training_images is:');
% disp(training_images);
% 
% 
% % Initialize the weights matrix
% weights_matrix_debug = zeros(num_pixels, num_pixels);
% 
% % Compute the weights matrix using the Hebbian Learning Rule
% for i = 1:num_pixels
%     for j = 1:num_pixels
%         if i ~= j
%             % Compute the weight between pixel i and pixel j
%             for col=1:3
%                 weights_matrix_debug(i, j) = weights_matrix_debug(i, j) + training_images(i, col) .* training_images(j, col);
%             end
%         end
%     end
% end
% 
% weights_matrix_debug = weights_matrix_debug/4;
% 
% % Display the resulting weights matrix
% disp('Weights matrix:');
% disp(weights_matrix_debug);
% 
% % Save the weights matrix to a file
% save('weights_matrix_debug.mat', 'weights_matrix_debug');
% 
% imagesc(weights_matrix_debug);
