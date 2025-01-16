%% Calculation Weights

% text files of images
file_names = {'number1.txt', 'number2.txt', 'number6.txt', 'number8.txt', 'number9.txt', 'number7.txt', 'number4.txt'};

% set matrix to zero
patterns = []; % 패턴을 저장할 배열

% read data from the text files
for f = 1:length(file_names)
    file_name = file_names{f};
    
    fprintf('filename is: [%s]\n', join(string(file_name), ', '));

    % read text file
    data = readmatrix(file_name);
    
    % save in 3D matrix format (row, column, number of patterns)
    patterns = cat(3, patterns, data);
    fprintf('data is: [%s]\n', join(string(data), ', '));
end


fprintf('patterns is: [%s] \n', join(string(patterns), ', '));


% data
[n_rows, n_cols, m_patterns] = size(patterns); % 행, 열, 패턴 수
n_pixels = n_rows * n_cols; % 총 픽셀 수

% matrix to zero (행과 열 기반)
C = zeros(n_rows, n_cols);

% Coupling Matrix calculatino 
for i = 1:n_rows
    for j = 1:n_cols
        % exclude self-calculation
        if i ~= j
            % Hebbian learning rule
            C(i, j) = sum(squeeze(patterns(i, j, :))) / n_pixels;
        end
    end
end

%% 
disp('Coupling Matrix C (row-column relationship):');
disp(C);

%% Finding Resistors Values

fprintf('Size of weights is: %d\n', size(C));

resistors = [];

for i = 1:n_rows
    for j = 1:n_cols
        if i ~= j
            resistors((i-1)*n_cols+j) = (C(i,j)*1000);
        end
    end 
end

%fprintf(['resistors are: ' repmat(' %0.000000001f ',1,numel(resistors)) '\n'],resistors);
fprintf('resistors are: [%s]\n', join(string(resistors), ', '));
fprintf('length of resistors is: %d\n', length(resistors));

%% 30x30 Weights Matrix

%clear all

% Define image dimensions and number of training images
rows = 5;
cols = 6;
num_images = 7;
num_pixels = rows * cols;

% Initialize a matrix to store the flattened training images
training_images = zeros(num_pixels, num_images);

filename = {'number1.txt', 'number2.txt', 'number6.txt', 'number8.txt', 'number9.txt', 'number7.txt', 'number4.txt'};


% Load 7 binary images
for i = 1:7
    % Replace 'imageX.png' with your binary image filenames
    %filename = sprintf('number%d.txt', k);
    %img = imread(filename{i}); % Read the binary image
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
            for col=1:7
                weights_matrix(i, j) = weights_matrix(i, j) + (training_images(i, col) .* training_images(j, col)) / num_pixels;
            end
            %weights_matrix(i, j) = sum(training_images(i, :) .* training_images(j, :)) / num_pixels;
        end
    end
end

% Display the resulting weights matrix
disp('Weights matrix:');
disp(weights_matrix);

% Optional: Save the weights matrix to a file
save('weights_matrix.mat', 'weights_matrix');


%% Finding the Resistors Values %%

% Extract the upper right triangle (excluding the diagonal)
upper_triangle_indices = triu(true(num_pixels), 1); % Create a logical mask for upper triangle
weights_vector = weights_matrix(upper_triangle_indices); % Extract the values

% Display the vector of weights
disp('Upper right triangle weights vector:');
disp(weights_vector);

% Optional: Save the weights matrix to a file
save('weights_matrix.mat', 'weights_matrix');
save('weights_vector.mat', 'weights_vector'); % Save the vector

%% Map to resistors %%

% Linear Mapping %%%%%%%%%%%%%%%%%%%%%
linear_res = weights_vector.*10000;

%disp(linear_res)
fprintf('Resistors Linear Mapping of Weights are: \n [%s]\n', join(string(linear_res), ', '));

% Optional: Save the Linear-transformed weights vector to a file
save('linear_res.mat', 'linear_res');



% Sigmoid Mappin %%%%%%%%%%%%%%%%%%%%
sig_res =  1 ./ (1 + exp(-weights_vector));

% Display the resulting sigmoid-transformed weights vector
fprintf('Resistors Sigmoid Mapping of Weights are: \n [%s]\n', join(string(sig_res), ', '));

% Optional: Save the sigmoid-transformed weights vector to a file
save('sigmoid_weights_vector.mat', 'sig_res');



% Exponential Mappin %%%%%%%%%%%%%%%%%%%%
exp_res =  exp(weights_vector.*10);

% Display the resulting sigmoid-transformed weights vector
fprintf('Resistors Exponential Mapping of Weights are: \n [%s]\n', join(string(exp_res), ', '));

% Optional: Save the sigmoid-transformed weights vector to a file
save('exp_weights_vector.mat', 'exp_res');



% Tuned Linear Mapping %%%%%%%%%%%%%%%%%%%%%
tune_linear_res = weights_vector.*100000 + 50000;

%disp(linear_res)
fprintf('Resistors Tuned Linear Mapping of Weights are: \n [%s]\n', join(string(tune_linear_res), ', '));

% Optional: Save the Linear-transformed weights vector to a file
save('tune_linear_res.mat', 'tune_linear_res');