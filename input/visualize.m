function images = visualize(filename)



    output = ReadData(filename);
    [n1, n2, n3, n4] = size(output);


    %y = ifftn(output) * sqrt(n1 * n2 * n3 * n4);
    y = ifftn(output) * sqrt(n3 * n4);

    images = uint8(absaa(y));


end