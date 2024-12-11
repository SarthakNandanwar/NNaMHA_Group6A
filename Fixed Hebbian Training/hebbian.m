{\rtf1\ansi\ansicpg1252\cocoartf2820
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;\red15\green112\blue16;\red0\green0\blue0;\red148\green0\blue242;
\red8\green0\blue255;}
{\*\expandedcolortbl;;\cssrgb\c0\c50196\c7451;\cssrgb\c0\c0\c0;\cssrgb\c65490\c3529\c96078;
\cssrgb\c5490\c0\c100000;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs20 \cf2 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 %% Calculation Weights\cf0 \strokec3 \
\
\cf2 \strokec2 % text files of images\cf0 \strokec3 \
file_names = \{\cf4 \strokec4 'number1.txt'\cf0 \strokec3 , \cf4 \strokec4 'number2.txt'\cf0 \strokec3 , \cf4 \strokec4 'number6.txt'\cf0 \strokec3 , \cf4 \strokec4 'number8.txt'\cf0 \strokec3 , \cf4 \strokec4 'number9.txt'\cf0 \strokec3 , \cf4 \strokec4 'number7.txt'\cf0 \strokec3 , \cf4 \strokec4 'number4.txt'\cf0 \strokec3 \};\
\
\cf2 \strokec2 % set matrix to zero\cf0 \strokec3 \
patterns = []; \cf2 \strokec2 % \uc0\u54056 \u53556 \u51012  \u51200 \u51109 \u54624  \u48176 \u50676 \cf0 \strokec3 \
\
\cf2 \strokec2 % read data from the text files\cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf5 \strokec5 for \cf0 \strokec3 f = 1:length(file_names)\
    file_name = file_names\{f\};\
    \
    fprintf(\cf4 \strokec4 'filename is: [%s]\\n'\cf0 \strokec3 , join(string(file_name), \cf4 \strokec4 ', '\cf0 \strokec3 ));\
\
    \cf2 \strokec2 % read text file\cf0 \strokec3 \
    data = readmatrix(file_name);\
    \
    \cf2 \strokec2 % save in 3D matrix format (row, column, number of patterns)\cf0 \strokec3 \
    patterns = cat(3, patterns, data);\
    fprintf(\cf4 \strokec4 'data is: [%s]\\n'\cf0 \strokec3 , join(string(data), \cf4 \strokec4 ', '\cf0 \strokec3 ));\
\cf5 \strokec5 end\cf0 \strokec3 \
\
\
fprintf(\cf4 \strokec4 'patterns is: [%s] \\n'\cf0 \strokec3 , join(string(patterns), \cf4 \strokec4 ', '\cf0 \strokec3 ));\
\
\
\pard\pardeftab720\partightenfactor0
\cf2 \strokec2 % data\cf0 \strokec3 \
[n_rows, n_cols, m_patterns] = size(patterns); \cf2 \strokec2 % \uc0\u54665 , \u50676 , \u54056 \u53556  \u49688 \cf0 \strokec3 \
n_pixels = n_rows * n_cols; \cf2 \strokec2 % \uc0\u52509  \u54589 \u49472  \u49688 \cf0 \strokec3 \
\
\cf2 \strokec2 % matrix to zero (\uc0\u54665 \u44284  \u50676  \u44592 \u48152 )\cf0 \strokec3 \
C = zeros(n_rows, n_cols);\
\
\cf2 \strokec2 % Coupling Matrix calculatino \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf5 \strokec5 for \cf0 \strokec3 i = 1:n_rows\
    \cf5 \strokec5 for \cf0 \strokec3 j = 1:n_cols\
        \cf2 \strokec2 % exclude self-calculation\cf0 \strokec3 \
        \cf5 \strokec5 if \cf0 \strokec3 i ~= j\
            \cf2 \strokec2 % Hebbian learning rule\cf0 \strokec3 \
            C(i, j) = sum(squeeze(patterns(i, j, :))) / n_pixels;\
        \cf5 \strokec5 end\cf0 \strokec3 \
    \cf5 \strokec5 end\cf0 \strokec3 \
\cf5 \strokec5 end\cf0 \strokec3 \
\
\pard\pardeftab720\partightenfactor0
\cf2 \strokec2 %% \cf0 \strokec3 \
disp(\cf4 \strokec4 'Coupling Matrix C (row-column relationship):'\cf0 \strokec3 );\
disp(C);\
\
\cf2 \strokec2 %% Finding Resistors Values\cf0 \strokec3 \
\
fprintf(\cf4 \strokec4 'Size of weights is: %d\\n'\cf0 \strokec3 , size(C));\
\
resistors = [];\
\
\pard\pardeftab720\partightenfactor0
\cf5 \strokec5 for \cf0 \strokec3 i = 1:n_rows\
    \cf5 \strokec5 for \cf0 \strokec3 j = 1:n_cols\
        \cf5 \strokec5 if \cf0 \strokec3 i ~= j\
            resistors((i-1)*n_cols+j) = (C(i,j)*1000);\
        \cf5 \strokec5 end\cf0 \strokec3 \
    \cf5 \strokec5 end \cf0 \strokec3 \
\cf5 \strokec5 end\cf0 \strokec3 \
\
\pard\pardeftab720\partightenfactor0
\cf2 \strokec2 %fprintf(['resistors are: ' repmat(' %0.000000001f ',1,numel(resistors)) '\\n'],resistors);\cf0 \strokec3 \
fprintf(\cf4 \strokec4 'resistors are: [%s]\\n'\cf0 \strokec3 , join(string(resistors), \cf4 \strokec4 ', '\cf0 \strokec3 ));\
fprintf(\cf4 \strokec4 'length of resistors is: %d\\n'\cf0 \strokec3 , length(resistors));\
\
\cf2 \strokec2 %% 30x30 Weights Matrix\cf0 \strokec3 \
\
\cf2 \strokec2 % Define image dimensions and number of training images\cf0 \strokec3 \
rows = 5;\
cols = 6;\
num_images = 7;\
num_pixels = rows * cols;\
\
\cf2 \strokec2 % Initialize a matrix to store the flattened training images\cf0 \strokec3 \
training_images = zeros(num_pixels, num_images);\
\
filename = \{\cf4 \strokec4 'number1.txt'\cf0 \strokec3 , \cf4 \strokec4 'number2.txt'\cf0 \strokec3 , \cf4 \strokec4 'number6.txt'\cf0 \strokec3 , \cf4 \strokec4 'number8.txt'\cf0 \strokec3 , \cf4 \strokec4 'number9.txt'\cf0 \strokec3 , \cf4 \strokec4 'number7.txt'\cf0 \strokec3 , \cf4 \strokec4 'number4.txt'\cf0 \strokec3 \};\
\
\
\cf2 \strokec2 % Load 7 binary images\cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf5 \strokec5 for \cf0 \strokec3 i = 1:7\
    \cf2 \strokec2 % Replace 'imageX.png' with your binary image filenames\cf0 \strokec3 \
    \cf2 \strokec2 %filename = sprintf('number%d.txt', k);\cf0 \strokec3 \
    \cf2 \strokec2 %img = imread(filename\{i\}); % Read the binary image\cf0 \strokec3 \
    data = readmatrix(filename\{i\});\
    disp(\cf4 \strokec4 'data is:'\cf0 \strokec3 );\
    disp(data);\
    img = imbinarize(data);  \cf2 \strokec2 % Ensure binary format (optional)\cf0 \strokec3 \
    img_vector = img(:);    \cf2 \strokec2 % Flatten the image into a vector\cf0 \strokec3 \
    training_images(:, i) = img_vector;\
\cf5 \strokec5 end\cf0 \strokec3 \
\
\pard\pardeftab720\partightenfactor0
\cf2 \strokec2 % Initialize the weights matrix\cf0 \strokec3 \
weights_matrix = zeros(num_pixels, num_pixels);\
\
\cf2 \strokec2 % Compute the weights matrix using the Hebbian Learning Rule\cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf5 \strokec5 for \cf0 \strokec3 i = 1:num_pixels\
    \cf5 \strokec5 for \cf0 \strokec3 j = 1:num_pixels\
        \cf5 \strokec5 if \cf0 \strokec3 i ~= j\
            \cf2 \strokec2 % Compute the weight between pixel i and pixel j\cf0 \strokec3 \
            \cf5 \strokec5 for \cf0 \strokec3 col=1:7\
                weights_matrix(i, j) = weights_matrix(i, j) + (training_images(i, col) .* training_images(j, col)) / num_pixels;\
            \cf5 \strokec5 end\cf0 \strokec3 \
            \cf2 \strokec2 %weights_matrix(i, j) = sum(training_images(i, :) .* training_images(j, :)) / num_pixels;\cf0 \strokec3 \
        \cf5 \strokec5 end\cf0 \strokec3 \
    \cf5 \strokec5 end\cf0 \strokec3 \
\cf5 \strokec5 end\cf0 \strokec3 \
\
\pard\pardeftab720\partightenfactor0
\cf2 \strokec2 % Display the resulting weights matrix\cf0 \strokec3 \
disp(\cf4 \strokec4 'Weights matrix:'\cf0 \strokec3 );\
disp(weights_matrix);\
\
\cf2 \strokec2 % Optional: Save the weights matrix to a file\cf0 \strokec3 \
save(\cf4 \strokec4 'weights_matrix.mat'\cf0 \strokec3 , \cf4 \strokec4 'weights_matrix'\cf0 \strokec3 );\
\
\
\cf2 \strokec2 %% Finding the Resistors Values %%\cf0 \strokec3 \
\
\cf2 \strokec2 % Extract the upper right triangle (excluding the diagonal)\cf0 \strokec3 \
upper_triangle_indices = triu(true(num_pixels), 1); \cf2 \strokec2 % Create a logical mask for upper triangle\cf0 \strokec3 \
weights_vector = weights_matrix(upper_triangle_indices); \cf2 \strokec2 % Extract the values\cf0 \strokec3 \
\
\cf2 \strokec2 % Display the vector of weights\cf0 \strokec3 \
disp(\cf4 \strokec4 'Upper right triangle weights vector:'\cf0 \strokec3 );\
disp(weights_vector);\
\
\cf2 \strokec2 % Optional: Save the weights matrix to a file\cf0 \strokec3 \
save(\cf4 \strokec4 'weights_matrix.mat'\cf0 \strokec3 , \cf4 \strokec4 'weights_matrix'\cf0 \strokec3 );\
save(\cf4 \strokec4 'weights_vector.mat'\cf0 \strokec3 , \cf4 \strokec4 'weights_vector'\cf0 \strokec3 ); \cf2 \strokec2 % Save the vector\cf0 \strokec3 \
\
\cf2 \strokec2 %% Map to resistors %%\cf0 \strokec3 \
\
\cf2 \strokec2 % Linear Mapping %%%%%%%%%%%%%%%%%%%%%\cf0 \strokec3 \
linear_res = weights_vector.*10000;\
\
\cf2 \strokec2 %disp(linear_res)\cf0 \strokec3 \
fprintf(\cf4 \strokec4 'Resistors Linear Mapping of Weights are: \\n [%s]\\n'\cf0 \strokec3 , join(string(linear_res), \cf4 \strokec4 ', '\cf0 \strokec3 ));\
\
\cf2 \strokec2 % Optional: Save the Linear-transformed weights vector to a file\cf0 \strokec3 \
save(\cf4 \strokec4 'linear_res.mat'\cf0 \strokec3 , \cf4 \strokec4 'linear_res'\cf0 \strokec3 );\
\
\
\
\cf2 \strokec2 % Sigmoid Mappin %%%%%%%%%%%%%%%%%%%%\cf0 \strokec3 \
sig_res =  1 ./ (1 + exp(-weights_vector));\
\
\cf2 \strokec2 % Display the resulting sigmoid-transformed weights vector\cf0 \strokec3 \
fprintf(\cf4 \strokec4 'Resistors Sigmoid Mapping of Weights are: \\n [%s]\\n'\cf0 \strokec3 , join(string(sig_res), \cf4 \strokec4 ', '\cf0 \strokec3 ));\
\
\cf2 \strokec2 % Optional: Save the sigmoid-transformed weights vector to a file\cf0 \strokec3 \
save(\cf4 \strokec4 'sigmoid_weights_vector.mat'\cf0 \strokec3 , \cf4 \strokec4 'sig_res'\cf0 \strokec3 );\
\
\
\
\cf2 \strokec2 % Exponential Mappin %%%%%%%%%%%%%%%%%%%%\cf0 \strokec3 \
exp_res =  exp(weights_vector.*10);\
\
\cf2 \strokec2 % Display the resulting sigmoid-transformed weights vector\cf0 \strokec3 \
fprintf(\cf4 \strokec4 'Resistors Exponential Mapping of Weights are: \\n [%s]\\n'\cf0 \strokec3 , join(string(exp_res), \cf4 \strokec4 ', '\cf0 \strokec3 ));\
\
\cf2 \strokec2 % Optional: Save the sigmoid-transformed weights vector to a file\cf0 \strokec3 \
save(\cf4 \strokec4 'exp_weights_vector.mat'\cf0 \strokec3 , \cf4 \strokec4 'exp_res'\cf0 \strokec3 );\
\
\
\
\cf2 \strokec2 % Tuned Linear Mapping %%%%%%%%%%%%%%%%%%%%%\cf0 \strokec3 \
tune_linear_res = weights_vector.*100000 + 50000;\
\
\cf2 \strokec2 %disp(linear_res)\cf0 \strokec3 \
fprintf(\cf4 \strokec4 'Resistors Tuned Linear Mapping of Weights are: \\n [%s]\\n'\cf0 \strokec3 , join(string(tune_linear_res), \cf4 \strokec4 ', '\cf0 \strokec3 ));\
\
\cf2 \strokec2 % Optional: Save the Linear-transformed weights vector to a file\cf0 \strokec3 \
save(\cf4 \strokec4 'tune_linear_res.mat'\cf0 \strokec3 , \cf4 \strokec4 'tune_linear_res'\cf0 \strokec3 );\
}