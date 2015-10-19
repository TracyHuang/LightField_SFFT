function transferData(datafile, outputfile, indices)

fid = fopen(datafile, 'r');

if fid < 0
    fprintf(2, 'Cannot open data file.\n');
end


images = importdata(datafile);
%data = fft2(imimages);

n_y = size(images, 1);
n_x = size(images, 2);
n_u = size(images, 3);
n_v = size(images, 4);
fprintf(1, '%d, %d, %d, %d\n', n_y, n_x, n_u, n_v);
data = zeros(size(images));

for i = 1 : n_y
    for j = 1 : n_x
        data(i, j, :, :) = fft2(squeeze(images(i, j, :, :))) / sqrt(n_u * n_v);
    end
end

%if(exist('indices', 'var'))
    
WriteData(outputfile, data);

%end



end
