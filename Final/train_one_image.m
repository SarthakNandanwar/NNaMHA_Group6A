%% Test %%%%%%%%%%%%%%%%%%%%%%%%% 

clear

% Define image dimensions and number of training images
rows_test = 7;
cols_test = 7;
num_images_test = 1;
num_pixels_test = rows_test * cols_test;

% Initialize a matrix to store the flattened training images
training_images_test = zeros(num_pixels_test, num_images_test);

filename_test = {'number0.txt'};


% Load 1 binary images
data_test = readmatrix(filename_test{1});
disp('data is:');
disp(data_test);
img_test = imbinarize(data_test);  % Ensure binary format (optional)
img_vector_test = img_test(:);    % Flatten the image into a vector
training_images_test(:, 1) = img_vector_test;

% Initialize the weights matrix
weights_matrix_test = zeros(num_pixels_test, num_pixels_test);

% Compute the weights matrix using the Hebbian Learning Rule
for i = 1:num_pixels_test
    for j = 1:num_pixels_test
        if i ~= j
            % Compute the weight between pixel i and pixel j
            
            weights_matrix_test(i, j) = weights_matrix_test(i, j) + training_images_test(i, 1) .* training_images_test(j, 1);
            
        end
    end
end

weights_matrix_test = weights_matrix_test/num_pixels_test;

% Display the resulting weights matrix
disp('Weights matrix:');
disp(weights_matrix_test);

% Save the weights matrix to a file
save('weights_matrix_test.mat', 'weights_matrix_test');

imagesc(weights_matrix_test);

%% Extracting the Weights Vector (Upper Triangle og the Weights Matrix) %%

% Extract the upper right triangle (excluding the diagonal)
upper_tri = triu(true(num_pixels_test), 1); % Create a logical mask for upper triangle
weights_vec = weights_matrix_test(upper_tri); % Extract the values

% Display the vector of weights
disp('Upper right triangle weights vector:');
disp(weights_vec);

% Save the weights matrix and vector to a file
% save('weights_matrix.mat', 'weights_matrix');
save('weights_vec.mat', 'weights_vec'); % Save the vector

%% Medium Linear Mapping %%%%%%%%%%%%%%%%%%%%%
med_linear_res_test = weights_vec.*6000000;
med_linear_res_mat_test = weights_matrix_test.*6000000 + 12000;
minValue_test = min(med_linear_res_test(:)); % Overall minimum
maxValue_test = max(med_linear_res_test(:)); % Overall maximum

for i=1:length(med_linear_res_test)
    if med_linear_res_test(i)==0
        med_linear_res_test(i) = maxValue_test;
    end
end

minValue_test2 = min(med_linear_res_test(:)); % Overall minimum

fprintf('The Min and Max of Medium Linear are: (%d, %d)\n', minValue_test2, maxValue_test);

% disp(linear_res)
% fprintf('Resistors Medium Linear Mapping of Weights are: \n [%s]\n', join(string(med_linear_res), ', '));
%fprintf('Resistors Matrix Medium Linear Mapping of Weights are: \n [%s]\n', join(string(med_linear_res_mat), ', '));

% Optional: Save the Linear-transformed weights vector to a file
% save('med_linear_res.mat', 'med_linear_res');
save('med_linear_res_mat_test.mat', 'med_linear_res_mat_test');
imagesc(med_linear_res_mat_test);