%This is script to crop cores using as input a CSV file with coordinates
%Remove the comments to also store as .mat files the coordinates

basePath = 'D:\users\fperez\NKI_TMAs_AF\';
outputfolder2='dearray\selected_channels_recroped_cords'; %To store the .mat files
outputfolder='dearray\selected_channels_recroped';
omePath = 'registration';
omeSuffix = '.ome.tif';
cropCoordsFileName = 'dearray\Recroped\Cores_to_recrop.csv';

%Selection of chanels
chanelSelected= [1,2,3,4,6,10,12,14,15,16,22,38,39,40,42,43,44,46,47,48];
numChannelNames = length(chanelSelected); % this is only to allocate memory before

%Select all samples
sampleList = dir( [ basePath 'TMA*' ] );

%Just for the next slides
list_of_Samples = [1 , 2 , 3, 4, 5, 8, 9, 10];

for sample = list_of_Samples
    sampleName = sampleList(sample).name;
    disp(sampleName);

    %sampleImage =  bfGetReader(strcat(basePath, sampleName, filesep, omePath, filesep, sampleName, omeSuffix));

    cropCoordstable = readtable(strcat(basePath, filesep, sampleName, filesep, cropCoordsFileName));
    [rows, columns] = size(cropCoordstable);
    coresFolder = strcat(basePath, filesep, sampleName, filesep, outputfolder);
    coordinatesFolder = strcat(basePath, filesep, sampleName, filesep, outputfolder2);
    mkdir(coresFolder);
    mkdir(coordinatesFolder);
    %For each of the cores in the CSV file crop them
    for x = 1:rows
        coreName = table2cell(cropCoordstable(x, 'Var1'));
        boundingBox = [cropCoordstable.Var5(x) + 1, cropCoordstable.Var6(x) + 1, cropCoordstable.Var7(x) + 1,  cropCoordstable.Var8(x) + 1];
        core = zeros(boundingBox(4), boundingBox(3), numChannelNames);
        rect = [cropCoordstable.Var5(x) + 1, cropCoordstable.Var6(x) + 1, cropCoordstable.Var7(x) + cropCoordstable.Var5(x) + 1,  cropCoordstable.Var8(x) + cropCoordstable.Var6(x) + 1];
        save(char(strcat(coordinatesFolder, filesep, replace(char(coreName), 'Core', ''), '_cropCoords.mat')), 'rect')
        
        box = num2cell(uint16(boundingBox));
        
        for iChan=1:length(chanelSelected)
            chan=chanelSelected(iChan);
            % Crop core from ome.tif
            core(:,:,iChan) = bfGetPlane(sampleImage, chan, box{:});
        end
        cc = uint16(core);
        try
            bfsave(cc,char(strcat(coresFolder, filesep, coreName, '.tif')), 'dimensionOrder', 'XYCZT');
        catch
            fprintf('Error detected');
        end
    end
end