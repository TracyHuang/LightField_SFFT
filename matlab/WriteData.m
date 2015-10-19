function WriteData(filename, data, indices, truncate_size)

    fid = fopen(filename, 'w');
    
    if fid < 0
        fprintf(2, 'Open %s failed.\n', filename);
    end
    
    if exist('indices', 'var')
        [index_length num_slices] = size(indices);
    else
        index_length = 0;
    end
    
    data_dims = size(data);
    
    if length(data_dims) == 2
        data_dims = [data_dims 1];
    end
    
    if ~exist('truncate_size', 'var')
        truncate_size = data_dims(1 : 2);
    else
        data_dims(1 : 2) = truncate_size;
    end
    
    
%     slice_size = data_dims(1) * data_dims(2);
    
    if ~exist('num_slices', 'var')
        num_slices = prod(data_dims(3 : end));
    end
    
    % the header
    fwrite(fid, num_slices, 'uint32');
    fwrite(fid, index_length, 'uint32');
    fwrite(fid, length(data_dims), 'uint32');
    fwrite(fid, data_dims, 'uint32');
    
    % the data
    for i = 1 : num_slices
        if (index_length ~= 0)
            fwrite(fid, indices(:, i), 'uint32');
        end
        fwrite(fid, real(data(1 : truncate_size(1), 1 : truncate_size(2), i)), 'double');
        fwrite(fid, imag(data(1 : truncate_size(1), 1 : truncate_size(2), i)), 'double');
    end
    
    fclose(fid);

end

