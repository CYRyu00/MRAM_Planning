function shape = AIMtoShape(M, AIM)

    shape = zeros(2*M+1,M); 
    shape(M+1,1) = 1; 
    rows=zeros(M,1);
    cols=zeros(M,1);
    rows(1)=M+1; cols(1)=1;
    
    if isempty(AIM)
        
        return;
    end
    
    direction = [];
    
    for i = 1: M
        % serach non-zero term in AIM(i,:), node i
        non_zero_indices_row = find(AIM(i, :) ~= 0);
        non_zero_values_row = AIM(i, non_zero_indices_row);
        
        for j = 1: length(non_zero_indices_row)
            col_index = non_zero_indices_row(j); % Column index of the current non-zero element
            non_zero_indices_col = find(AIM(:, col_index) ~= 0);
            non_zero_values_col = AIM(non_zero_indices_col, col_index);
            %disp(non_zero_indices_col)
            for k = 1:length(non_zero_indices_col)
                row_idx = non_zero_indices_col(k);
                %disp(row_idx)
                if row_idx ~= i
                    a = AIM(row_idx,col_index);
                    if a == -1
                        direction = [1, 0];         
                    elseif a == +1
                        direction = [-1,0];
                    elseif a == -2
                        direction = [0,1];
                    elseif a == +2
                        direction = [0,-1];
                    end
                    pos = [rows(i), cols(i)] + direction;
                    rows(row_idx) = pos(1); cols(row_idx) = pos(2);
                    shape(rows(row_idx) , cols(row_idx) ) = row_idx;
                end           
            end
        end
    end
%TODO labeled - > 1 1 1 1
end
