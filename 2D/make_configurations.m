function result = make_configurations(num_AMs, max_L)
    % This is the main function that calls the helper function findCombinations
    result = findCombinations(num_AMs, max_L);
end

function valid_L = findCombinations(num_AMs, max_L)
    % Helper function to find combinations of partitions
    all_combinations = partition(num_AMs, max_L);
    valid_L = {};
    
    % Filter combinations to meet the conditions
    for i = 1:length(all_combinations)
        L = all_combinations{i};
        if sum(L) == num_AMs && all(L > 0) && length(L) <= max_L
            valid_L{end+1} = L; %#ok<AGROW> Append valid combination
        end
    end
end

function partitions = partition(n, max_parts)
    % Recursive function to find all partitions of n into up to max_parts parts
    if n == 0
        partitions = {[]};
    elseif max_parts == 0 || n < 0
        partitions = {};
    else
        partitions = {};
        for k = 1:n
            sub_partitions = partition(n - k, max_parts - 1);
            for j = 1:length(sub_partitions)
                if length([k, sub_partitions{j}]) <= max_parts  % Limit length directly
                    partitions{end+1} = [k, sub_partitions{j}]; %#ok<AGROW>
                end
            end
        end
    end
end
