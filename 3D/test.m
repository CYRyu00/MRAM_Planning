m =4;
result = generateAllShapes(m);
%disp(result(:,:,:)); % Displays the size of the 3D array

%AIM_Test = generateAIM(result(:,:,5));
result_2 = labelNodes(result(:,:,7));
disp(result_2);
max(result_2,[],"all");

result_3 = generateAIM(result_2);
disp(result_3)

