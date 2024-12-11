%% Read Matrix

% Read the file as text
filename = 'number1.txt'; % Replace with the correct path
fileID = fopen(filename, 'r');
lines = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);

lines = lines{1}; % Extract the cell array of lines

% Initialize matrix
mat1 = zeros(5, 6); % Preallocate a 5x6 matrix
rowIdx = 0; % To track the current row

for i = 1:length(lines)
    line = strtrim(lines{i}); % Remove leading/trailing whitespace
    
    % Skip non-numeric lines
    if isempty(line) || startsWith(line, 'number')
        continue;
    end
    
    % Convert line to a numeric row
    row = str2num(line); %#ok<ST2NM>
    if length(row) ~= 6
        error('Row "%s" is not 6 elements long.', line);
    end
    
    % Add the row to the matrix
    rowIdx = rowIdx + 1;
    mat1(rowIdx, :) = row;
    
    % Stop after 5 rows (since it's a 5x6 matrix)
    if rowIdx == 5
        break;
    end
end

% Display the matrix
disp('Matrix for number1:');
disp(mat1);

weights = zeros(5, 6);
for i = 1:5
    for j = 1:6
        if i == j
            weights(i, j) = 0; % Set diagonal elements to zero
        else
            weights(i, j) = (1/30) * mat1(i, i) * mat1(i, j);
        end
    end
end

disp(weights);
