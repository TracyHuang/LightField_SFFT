function splitData(data, output_file_prefix, num_splits)

%fid = fopen(datafile, 'r');

%if fid < 0
%    fprintf(2, 'Cannot open data file.\n');
%end

%data = ReadData(datafile);

%images = importdata(datafile);
%data = fft2(imimages);

%n_x = size(images, 1);
%n_y = size(images, 2);
%n_u = size(images, 3);
%n_v = size(images, 4);
%data = zeros(size(images));

%for i = 1 : n_x
%    for j = 1 : n_y
%        data(i, j, :, :) = fft2(squeeze(images(i, j, :, :))) / sqrt(n_u * n_v);
%    end
%end

data_dims = size(data);
num_slices = prod(data_dims(3 : end));

%if mod(num_slices, num_splits) ~= 0
%   fprintf(2, 'error num_slices can not be divided by num_splits\n');
%   return
%end

num_slices_per_split = ceil(num_slices / num_splits);

for i = 1 : (num_splits - 1)
    if i <= 10
    output_file_name = strcat(output_file_prefix, '0',  num2str(i - 1), '.dat');
    else
    output_file_name = strcat(output_file_prefix, num2str(i - 1), '.dat');    
    end
    fid = fopen(output_file_name, 'w');
    
    if fid < 0
       fprintf(2, 'Open %s failed.\n', output_file_name);
    end
    
    fwrite(fid, num_slices_per_split, 'uint32');
    fwrite(fid, 0, 'uint32');
    fwrite(fid, length(data_dims) - 1, 'uint32');
    fwrite(fid, [size(data, 1), size(data, 2), num_slices_per_split], 'uint32');
    
    for j = ((i-1) * num_slices_per_split + 1) : i * num_slices_per_split
       fwrite(fid, real(data(1 : data_dims(1), 1 : data_dims(2), j)), 'double');
       fwrite(fid, imag(data(1 : data_dims(1), 1 : data_dims(2), j)), 'double'); 
    end
    fclose(fid);
end

% special case the last one
last_index = num_splits - 1;
if last_index <= 10
    output_file_name = strcat(output_file_prefix, '0', num2str(last_index), '.dat');
else
    output_file_name = strcat(output_file_prefix, num2str(last_index), '.dat');
end
    
fid = fopen(output_file_name, 'w');

if fid < 0
   fprintf(2, 'Open %s failed\n', output_file_name);
end

num_slices_left = num_slices - num_slices_per_split * (num_splits - 1);

fwrite(fid, num_slices_left, 'uint32');
fwrite(fid, 0, 'uint32');
fwrite(fid, length(data_dims) - 1, 'uint32');
fwrite(fid, [size(data, 1), size(data, 2), num_slices_left], 'uint32');

for j = (num_slices_per_split * (num_splits - 1) + 1) : num_slices
    fwrite(fid, real(data(1 : data_dims(1), 1 : data_dims(2), j)), 'double');
    fwrite(fid, imag(data(1 : data_dims(1), 1 : data_dims(2), j)), 'double'); 
end
    
end