function output = catImages(image1, image2, n1, n2, n3, n4)

    iwidth = max(size(image1,2),size(image2,2));
    if size(image1,2) < iwidth
        image1(1,iwidth,1) = 0;
    end
    if size(image2,2) < iwidth
        image2(1,iwidth,1) = 0;
    end
    output = cat(1,image1,image2);



end