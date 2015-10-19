function output_movie = result2movie(dir1, dir2, dir3,  file_name)

    %m = colormap();
    fileList1 = getAllFiles(dir1);
    fileList2 = getAllFiles(dir2);
    fileList3 = getAllFiles(dir3);
    t=1;
    for i = 1 : 17
       for j = 1 : 17
       filename = fileList1(t);
       image1 = imread(filename{:});
       filename = fileList2(t);
       image2 = imread(filename{:});
       filename = fileList3(t);
       image3 = imread(filename{:});
       image = con3pic(image1, image2, image3);
       %F(t) = im2frame(squeeze(images(i, j, :, :)), gray(256));
       F(t) = im2frame(image);
       t = t + 1;
       end
        
    end
    
    
     output_movie = F;

    movie2avi(F, file_name);

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


function image = con3pic(image1, image2, image3)

iwidth = max(size(image1,2),size(image2,2));
iwidth = max(iwidth, size(image3,2));
if size(image1,2) < iwidth
  image1(1,iwidth,1) = 0;
end
if size(image2,2) < iwidth
  image2(1,iwidth,1) = 0;
end
if size(image3, 2) < iwidth
   image3(1,iwidth,1) = 0; 
end
image = cat(2,image1,image2, image3);



end