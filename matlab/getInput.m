function getInput(directory, n1, n2, n3, n4, output_file_prefix, num_splits)

%this function is just a generanization of the process of generating input
% for details please check splitData.m and generateInput.m

input = generateInput(directory, n1, n2, n3, n4);

splitData(input, output_file_prefix, num_splits);

end