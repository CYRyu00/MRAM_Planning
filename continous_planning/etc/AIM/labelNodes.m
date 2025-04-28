function labeledMap = labelNodes(occupancyMap)
    % Input:
    %   occupancyMap - A 2D array representing the occupancy map (binary values)
    % Output:
    %   labeledMap - A 2D array where the 1s are labeled starting from the center (n+1, 1)

    % Get the size of the input map
    [rows, cols] = size(occupancyMap);
    
    % Ensure the grid dimensions match the expected structure
    if mod(rows, 2) == 0
        error('The number of rows must be odd for proper centering.');
    end
    
    % Determine the center point (n+1, 1)
    centerRow = (rows + 1) / 2;
    centerCol = 1;
    
    % Initialize the labeled map with zeros
    labeledMap = zeros(rows, cols);
    
    % Initialize the label counter
    label = 1;
    
    % Ensure the center point is labeled first
    if occupancyMap(centerRow, centerCol) == 1
        labeledMap(centerRow, centerCol) = label;
        label = label + 1;
    else
        error('The center point must be occupied (1).');
    end
    
    % Initialize a queue for BFS
    queue = [centerRow, centerCol];
    
    % Define movement directions: [Right, Down, Left, Up]
    directions = [0, 1; 1, 0; 0, -1; -1, 0];
    
    % Perform BFS to label all connected components
    while ~isempty(queue)
        % Dequeue the first element
        [currentRow, currentCol] = deal(queue(1, 1), queue(1, 2));
        queue(1, :) = [];
        
        % Explore all neighbors
        for d = 1:size(directions, 1)
            newRow = currentRow + directions(d, 1);
            newCol = currentCol + directions(d, 2);
            
            % Check bounds and if the cell is unvisited and occupied
            if newRow >= 1 && newRow <= rows && ...
               newCol >= 1 && newCol <= cols && ...
               occupancyMap(newRow, newCol) == 1 && ...
               labeledMap(newRow, newCol) == 0
           
                % Label the cell and enqueue it
                labeledMap(newRow, newCol) = label;
                label = label + 1;
                queue = [queue; newRow, newCol];
            end
        end
    end
end