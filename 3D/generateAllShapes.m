function occupancyCombinations = generateAllShapes(n)
    % Input:
    %   n - Determines the grid size [(2n+1) x n] and maximum number of 1s
    % Output:
    %   occupancyCombinations - 3D array with all valid 2D occupancy grid combinations

    % Define grid size
    rows = 2 * n + 1; % Height of the grid
    cols = n;         % Width of the grid

    % Initialize result as a cell array to handle variable size grids
    occupancyCombinations = {};
    
    % Recursive function to generate valid grids
    function generateGrids(grid, currentRow, currentCol, onesCount)
        % Stop recursion if the number of 1s exceeds n
        if onesCount > n
            return;
        end
        
        % Store the grid if it's valid (add to the combinations list)
        occupancyCombinations{end+1} = grid;
        
        % Explore all neighbors to add new 1s (connectivity constraint)
        directions = [0, 1; 1, 0; 0, -1; -1, 0]; % Right, Down, Left, Up
        for i = 1:size(directions, 1)
            newRow = currentRow + directions(i, 1);
            newCol = currentCol + directions(i, 2);
            
            % Ensure the new position is within bounds and unoccupied
            if newRow >= 1 && newRow <= rows && ...
               newCol >= 1 && newCol <= cols && ...
               grid(newRow, newCol) == 0
                % Create a new grid with the updated position
                newGrid = grid;
                newGrid(newRow, newCol) = 1;
                % Recur with the new grid
                generateGrids(newGrid, newRow, newCol, onesCount + 1);
            end
        end
    end

    % Start with an empty grid
    initialGrid = zeros(rows, cols);
    initialGrid(n+1, 1) = 1; % Ensure (n+1, 1) is always 1
    generateGrids(initialGrid, n+1, 1, 1); % Start recursive generation

    % Convert cell array to 3D array
    numGrids = numel(occupancyCombinations);
    occupancyCombinations3D = zeros(rows, cols, numGrids);
    for i = 1:numGrids
        occupancyCombinations3D(:, :, i) = occupancyCombinations{i};
    end
    % Output the result
    occupancyCombinations = occupancyCombinations3D;
end
