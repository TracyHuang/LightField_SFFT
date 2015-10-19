function [data, used] = mergeResult(directory, n1, n2, n3, n4)

fileList = getAllFiles(directory);
data = zeros(n1, n2, n3, n4);
slice_index = 0; % index of the current slice
used = zeros(n3, n4);
    for i = 1 : length(fileList)   
    file_name = fileList(i);
    %file_path = strcat(directory, file_name);
    file_data = ReadData(file_name{:});
    file_num_slices = size(file_data, 3);
        for j = 1 : file_num_slices
            slice_index = slice_index + 1;
            [index_3d, index_4d] = findIndex(slice_index, n3);
            fprintf(1, 'index is %d, %d\n', index_3d, index_4d);
            used(index_3d, index_4d) = used(index_3d, index_4d) + 1;
            for k = 1 : n1
                for l = 1 : n2
                    data(k, l, index_3d, index_4d) = file_data(k, l, j);
                    
                end
            end
        end
    end
    
    
end


function [index_3d, index_4d] = findIndex(slice_index, num_3d)

%index_3d = floor((slice_index - 1) / n4) + 1;           %row
%index_4d = slice_index - (index_3d - 1) * n4;    %column

index_4d = floor((slice_index - 1) / num_3d) + 1;
index_3d = slice_index -  (index_4d - 1) * num_3d;


end


function fileList = getAllFiles(dirName)

  dirData = dir(dirName);      %# Get the data for the current directory
  dirIndex = [dirData.isdir];  %# Find the index for directories
  fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
  if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                       fileList,'UniformOutput',false);
  end
  subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
                                               %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    fileList = [fileList; getAllFiles(nextDir)];  %# Recursively call getAllFiles
  end

end
