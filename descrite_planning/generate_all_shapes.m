function shapes = generate_all_shapes(m,K,L,core)
    shapes = cell(1,m);
    
    shape = zeros(K,L);
    shape(core(1),core(2)) = 2;
    shapes{1} = {shape};
    
    % Get size of shape
    [rows, cols] = size(shape);
    %labled_shape = zeros(rows, cols);
    
    % Define 4-directional neighbors (up, down, left, right)
    directions = [-1 0;  % Up
                   1 0;  % Down
                   0 -1; % Left
                   0 1]; % Right
    
    for i=1:m-1
        fprintf("Generate shape with %d AM\n", i+1);
        unique_shapes = {};
        for j=1:length(shapes{i})
            shape = shapes{i}{j};
            [one_rows, one_cols] = find(shape ~= 0);
            
            for k = 1:length(one_rows)
                r = one_rows(k);
                c = one_cols(k);
                for d = 1:size(directions, 1)
                    nr = r + directions(d, 1);
                    nc = c + directions(d, 2);
                    
                    % Check if neighbor is within shape bounds and not already 1
                    if nr >= 1 && nr <= rows && nc >= 1 && nc <= cols && shape(nr, nc) == 0
                        new_shape = shape; new_shape(nr, nc) = 1;
                        shapes{i+1}{end+1} = new_shape; % Add new_shape to shapes{i+1}
                    end
                end
            end
        end
        % Remove duplicate shapes
        shapes{i+1} = unique(cellfun(@(x) mat2str(x), shapes{i+1}, 'UniformOutput', false));
        shapes{i+1} = cellfun(@(x) eval(x), shapes{i+1}, 'UniformOutput', false);
    end

end
