% text files of images
file_names = {'number1.txt', 'number2.txt', 'number6.txt', 'number8.txt', 'number9.txt', 'number7.txt', 'number4.txt'};

% set matrix to zero
patterns = []; % 패턴을 저장할 배열

% read data from the text files
for f = 1:length(file_names)
    file_name = file_names{f};
    
    % read text file
    data = readmatrix(file_name);
    
    % save in 3D matrix format (row, column, number of patterns)
    patterns = cat(3, patterns, data);
end

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

% 
disp('Coupling Matrix C (row-column relationship):');
disp(C);
