function AIM = generateAIM(labeledMap)
    m=max(labeledMap,[],"all");
    uniqueLabels = 1:1:m;
    
    AIM=[];

    % Define movement directions: [Right, Down, Left, Up]
    directions = [0, 1; 1, 0; 0, -1; -1, 0];
    
    % Iterate over each label
    for label = uniqueLabels
        % Find the position of the current label
        [row, col] = find(labeledMap == label);
        
        % Initialize neighbor set
        neighborSet = [];
        
        % Check neighbors
        for d = 1:size(directions, 1)
            newRow = row + directions(d, 1);
            newCol = col + directions(d, 2);
            
            % Ensure the neighbor is within bounds
            if newRow >= 1 && newRow <= size(labeledMap, 1) && ...
               newCol >= 1 && newCol <= size(labeledMap, 2)
                neighborLabel = labeledMap(newRow, newCol);
                
                % Add to neighbors if it's non-zero and not the current label
                if neighborLabel ~= 0 && neighborLabel ~= label
                    neighborSet = unique([neighborSet, neighborLabel]);
                    
                    newAIM = zeros(m,1);
                    if directions(d,:) == [1, 0]
                        a = +1;
                    elseif directions(d,:) ==[-1,0]
                        a = -1;
                    elseif directions(d,:) == [0,1]
                        a = 2;
                    elseif directions(d,:) == [0,-1]
                        a = -2;
                    end
                    newAIM(label,1) = a;
                    newAIM(neighborLabel,1) = -a;

                    % Ensure `newAIM` has the same number of rows as `AIM`
                    if isempty(AIM) % If AIM is empty, initialize with `newAIM`
                        AIM = newAIM;
                    else
                        % Make sure `newAIM` is initialized completely before comparison
                        newAIM(label, 1) = a;
                        newAIM(neighborLabel, 1) = -a;
                    
                        % Check if `newAIM` is a duplicate
                        isDuplicate = any(all(AIM == newAIM, 1));
                    
                        % Add `newAIM` as a new column only if it is unique
                        if ~isDuplicate
                            AIM = [AIM, newAIM];
                        end
                    end

                    %AIM = [AIM ,newAIM];
                end
            end
        end
        
       
    end
end
