targetMatrix = ...
[   0	0	0	0	0	0
    0	0	0	0	0	0
    0	0	0	0	0	0
    0	0	0	0	0	0
    0	0	0	0	0	0
    0	0	0	0	0	0
    2	0	0	0	0	0
    1	0	0	0	0	0
    1	0	0	0	0	0
    1	0	0	0	0	0
    1	0	0	0	0	0
    0	0	0	0	0	0
    0	0	0	0	0	0];

foundIndex = [];
for i = 1:numel(shapes)
    if isequal(shapes{i}, targetMatrix)
        foundIndex = i;
        fprintf("%d, J*: %f\n", i, all_optimal_value{i})
        break;
    end
end