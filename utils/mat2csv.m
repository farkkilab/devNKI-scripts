basePath = 'D:\users\fperez\NKI_TMAs_AF\';
outputfolder='dearray\cropCoords\'; %To store the .csv files
cropCoordsPath = 'dearray\cropCoords\';
cropCoordsFileName = '*_cropCoords.mat';
csvCoordsFileSuffix = '_cropCoords.csv';


%Select all samples
sampleList = dir( [ basePath 'TMA*' ] );

%Just for the next samples
%list_of_Samples = [4, 5, 8, 10];
list_of_Samples = [1 : 10];

for sample = list_of_Samples
    sampleName = sampleList(sample).name;
    disp(sampleName);
    
    cordsdir = strcat(basePath, filesep, sampleName, filesep, cropCoordsPath, filesep);

    cropCoordsFiles = dir( [ basePath filesep sampleName filesep cropCoordsPath filesep cropCoordsFileName ] );

    for coreCoords = 1:length(cropCoordsFiles)
        coreCoordsName = cropCoordsFiles(coreCoords).name;
        splitName = strsplit(coreCoordsName, '_');
        iCore = splitName{1};
        fprintf('Sample: %s - core %s Started \n', sampleName, iCore);
        %Coordinate .mat files must contain a 'rect' object
        croppingdata = load( [ cropCoordsFiles(coreCoords).folder filesep coreCoordsName ] );
        rect = croppingdata.rect;
        writematrix(rect,char(strcat(cordsdir, iCore, csvCoordsFileSuffix)),'Delimiter',',');
    end
end