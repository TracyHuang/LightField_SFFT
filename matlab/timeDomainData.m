function data = timeDomainData(directory, n1, n2, n3, n4)

fileList = getAllFiles(directory);
data = zeros(n1, n2, n3, n4);

for i = 1 : length(fileList)
    file_name = fileList(i);
    image_file = imread(file_name{:});
    new_image = imresize(image_file, 0.5);
    %n3 = size(new_image, 3);
    %n4 = size(new_image, 4);
    
    if (size(new_image, 1) ~= n3 || size(new_image, 2) ~= n4)
        fprintf(2, 'error new image dimension dimensions do not match with argument\n');
        return
    end
    
    if length(size(new_image)) ~= 3
       fprintf(2, 'error image does not have 3 dimensions\n');
       return
    end
    grey_image = zeros(size(new_image, 1), size(new_image, 2));
    for j = 1 : size(grey_image, 1)
        for k = 1 : size(grey_image, 2)
            grey_image(j, k) = 0.257 * new_image(j, k, 1) + 0.504 * new_image(j, k, 2) + 0.098 * new_image(j, k, 3) + 16.0;
        end
    end
    
    row = floor((i - 1) / n2) + 1;
    column = i - (row - 1) * n2;
    
    for l = 1 : n3
        for m = 1 : n4
            data(row, column, l, m) = grey_image(l, m);
        end 
    end
    
end

end