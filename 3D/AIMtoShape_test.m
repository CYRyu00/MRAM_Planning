m = 5;
result = generateAllShapes(m);
%disp(result(:,:,:)); % Displays the size of the 3D array

%AIM_Test = generateAIM(result(:,:,5));
for i = 1:length(result())
    labeled = labelNodes(result(:,:,i));
    %disp(labeled);
    AIM = generateAIM(m, labeled);
    
    shape = AIMtoShape(m,AIM);
    
    if ~isequal(shape ,labeled)
        disp(isequal(shape ,labeled))
        disp(i)
    end
end

