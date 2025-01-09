
close all
clear


% Step 1: Define the binary matrix (values 0 and 1)
binaryMatrix = [
    1 0 1 0 1;
    0 1 0 1 0;
    1 0 1 0 1;
    0 1 0 1 0
];

% Step 2: Convert to logical format (optional, for clarity)
binaryImage = logical(binaryMatrix);

% Step 3: Display the binary image
imshow(binaryImage, 'InitialMagnification', 500); % Scaled up for visibility
title('Binary Image: Black and White');

% Step 4: Save the image as a binary PNG
imwrite(binaryImage, 'binaryImage.png');

disp('Image saved as binaryImage.png');