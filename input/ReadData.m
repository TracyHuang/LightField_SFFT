function [data, indices, effective_num] = ReadData(filename, subsampling)

    fid = fopen(filename, 'r');
   
    if fid < 0
        fprintf(2, 'Open %s failed.\n', filename);
        data = [];
        effective_num = 0;
        indices = [];
        return;
    end
    
    
    % first let's read the header
    num_slices = fread(fid, 1, 'uint32');
    
    if (isempty(num_slices))
        effective_num = 0;
        data = [];
        indices = [];
        return;
    end
    
    index_length = fread(fid, 1, 'uint32');
    num_data_dims = fread(fid, 1, 'uint32');
    data_dims = fread(fid, [1 num_data_dims], 'uint32');
    
    
    % read the data
    if ~exist('subsampling', 'var')
        data = zeros(data_dims);
    else
        data = zeros([size(subsampling, 2), data_dims(3 : end)]);
    end
    indices = zeros(index_length, num_slices);
    
    for i = 1 : num_slices
        if index_length ~= 0
            tmp_indices = fread(fid, [index_length 1], 'uint32');
            
            if length(tmp_indices) == index_length
                indices(:, i) = tmp_indices;
            else
                effective_num = i - 1;
                return;
            end
        end
        
        real_d = fread(fid, [data_dims(1) data_dims(2)], 'double');
        image_d = fread(fid, [data_dims(1) data_dims(2)], 'double');
        
        
        if (size(real_d, 1) == data_dims(1)) && (size(real_d, 2) == data_dims(2)) && (size(image_d, 1) == data_dims(1)) && (size(image_d, 2) == data_dims(2))
            if ~exist('subsampling', 'var')
                data(:, :, i) = real_d + 1i * image_d;
            else
                for subsample_index = 1 : size(subsampling)
                    data(subsample_index, i) = real_d(subsampling(1, subsample_index), subsampling(2, subsample_index)) + 1i * image_d(subsampling(1, subsample_index), subsampling(2, subsample_index));
                end
            end
        else
            effective_num = i - 1;
            return;
        end
    end
    
    effective_num = num_slices;
        
    fclose(fid);

end