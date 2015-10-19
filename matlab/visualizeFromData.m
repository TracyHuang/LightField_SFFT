function images = visualizeFromData(data)

    [n1, n2, n3, n4] = size(data);


    y = ifftn(data) * sqrt(n1 * n2 * n3 * n4);
    %y = ifftn(data) * sqrt(n3 * n4) * 14.7128;

    images = uint8(real(y));


end