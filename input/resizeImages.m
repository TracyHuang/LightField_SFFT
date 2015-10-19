function resizeImages(directory, n2) %n2 is number of columns

fileList = getAllFiles(directory);

for i = 1 : length(fileList)
   file_name = fileList(i);
   image_file = imread(file_name{:});
   new_image = imresize(image_file, 0.5);
   row = floor((i - 1) / n2);
   column = i - row * n2 - 1;
   if row < 10
       if column < 10
       new_image_name = strcat('resized_', '0', num2str(uint8(row)), '_0', num2str(uint8(column)), '.png');
       else
        new_image_name = strcat('resized_', '0', num2str(uint8(row)), '_', num2str(uint8(column)), '.png'); 
       end
   
   else
       if column < 10
       new_image_name = strcat('resized_', num2str(uint8(row)), '_0', num2str(uint8(column)), '.png');
       else
       new_image_name = strcat('resized_', num2str(uint8(row)), '_', num2str(uint8(column)), '.png');    
       end
   end
   imwrite(new_image, new_image_name);
end


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