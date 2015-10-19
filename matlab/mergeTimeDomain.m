function image = mergeTimeDomain(directory, n1, n2, n3, n4)

image = zeros(n1, n2, n3, n4);
fileList = getAllFiles(directory);

for i = 1 : length(fileList)   
    file_name = fileList(i);
    %file_path = strcat(directory, file_name);
    file_data = imread(file_name{:});
    
    index_1d = floor((i-1) / n1) + 1;
    index_2d = i - n2 * (index_1d - 1);
    for j = 1 : n3
        for k = 1 : n4
            image(index_1d, index_2d, j, k) = file_data(j, k);
        end
    end
    
end

image = uint8(image);

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