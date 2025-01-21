% Transform file.txt pattern into init capacitor voltage

% Define image dimensions and number of training images
rows = 7;
cols = 7;
num_images = 7;
num_pixels = rows * cols;

% Initialize a matrix to store the flattened training images
training_images = zeros(num_pixels, num_images);

%filename = {'number1.txt', 'number2.txt', 'number6.txt', 'number8.txt', 'number9.txt', 'number7.txt', 'number4.txt'};
% Get the current working directory
folder_path = pwd;

% Get a list of all .txt files in the current folder
files = dir(fullfile(folder_path, '*.txt'));

% Extract the names of the files
file_name = {files.name};

% Count the number of .txt files
num_files = length(file_name);

% Display the number of .txt files
fprintf('There are %d .txt files in the current folder.\n', num_files);

% Display the result
disp('File names:');
disp(file_name);

% Load 7 binary images
for i = 1:num_files

    data = readmatrix(file_name{i});
    disp('data is:');
    disp(data);
    img = imbinarize(data);  % Ensure binary format (optional)
    img_vector = img(:);    % Flatten the image into a vector
    training_images(:, i) = img_vector;
end
target_value_0 = 0.2;
target_value_1 = 0.6;
% Replace values using logical indexing
training_images(training_images == 0) = target_value_0;
training_images(training_images == 1) = target_value_1;
V_c_init_col_mat = training_images;
V_c_init_col_mat (1, :) = 0; 


save('V_c_init_col_mat.mat',"V_c_init_col_mat");
