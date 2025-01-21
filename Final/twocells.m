%% 2x1 Weights Training Matrix Calculation 

clear

% Define image dimensions and number of training images
rows = 1;
cols = 2;
num_images = 1;
num_pixels = rows * cols;

% Initialize a matrix to store the flattened training images
training_images = zeros(num_pixels, num_images);

filename = {'2cells.txt'};


data = readmatrix(filename{1});
disp('data is:');
disp(data);
img = imbinarize(data);  % Ensure binary format (optional)
img_vector = img(:);    % Flatten the image into a vector
training_images(:, 1) = img_vector;

% Initialize the weights matrix
weights_mat = zeros(num_pixels, num_pixels);

% Compute the weights matrix using the Hebbian Learning Rule
for i = 1:num_pixels
    for j = 1:num_pixels
        if i ~= j
            % Compute the weight between pixel i and pixel j
            
            weights_mat(i, j) = weights_mat(i, j) + training_images(i, 1) .* training_images(j, 1);
            
        end
    end
end

weights_mat = weights_mat/num_pixels;

% Display the resulting weights matrix
disp('Weights matrix:');
disp(weights_mat);

% Save the weights matrix to a file
save('2cells_wmat.mat', 'weights_mat');

imagesc(weights_mat);


%% Map to Resistors %%

% Linear Mapping %%%%%%%%%%%%%%%%%%%%%
twocells_rest_mat = weights_mat + 40000;


% Save the Linear-transformed weights matrix to a file
save('twocells_res_mat.mat', 'twocells_rest_mat');
imagesc(twocells_rest_mat);
 