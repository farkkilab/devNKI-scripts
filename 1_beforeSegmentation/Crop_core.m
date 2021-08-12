% Script used to extract relevant channels for labeling in ilastik and get
% the pixel probability maps
% Made by Fernando

basePath = 'D:\users\fperez\NKI_TMAs_AF\';
destPath='D:\users\fperez\NKI_TMAs_AF';
outputfolder='selected_Channels';
omePath = 'registration';
omeSuffix = '.ome.tif';
cropCoordsPath = 'dearray';
cropCoordsFileName = '*_cropCoords.mat';

%Chanel information
channelNames = readtable( [basePath filesep 'channel_list.csv'], 'ReadVariableNames', false);

%Selection of chanels
chanelSelected= [1,2,3,4,6,10,12,14,15,16,22,38,39,40,42,43,44,46,47,48];

%Extra DNA channels
%DNAs=5:4:41;

%Remove the extra DNA channels to reduce space in hard disk
%chanelSelected(DNAs)= [];

numChannelNames = length(chanelSelected); % this is only to allocate memory before loops


%chanelSelected = 1;

%Select random sample
sampleList = dir( [ basePath 'TMA*' ] );


for sample = 1:length(sampleList)
%for sample = 1
        %Number of crop files
        sampleName = sampleList(sample).name;
        disp(sampleName)
        tic
        %Read ometif image for that sample
        sampleImage =  bfGetReader( [ basePath sampleName filesep omePath filesep sampleName omeSuffix ] );
        numChannels = sampleImage.getImageCount();
        cropCoordsFiles = dir( [ basePath filesep sampleName filesep cropCoordsPath filesep cropCoordsFileName ] );


        for coreCoords = 1:length(cropCoordsFiles)

            %Get cores, crop them from big TIFF
            cropCoordsFiles = dir( [ basePath filesep sampleName filesep cropCoordsPath filesep cropCoordsFileName ] );

            coresFolder = [ destPath filesep sampleName filesep cropCoordsPath filesep outputfolder];
            mkdir(coresFolder);

            coreCoordsName = cropCoordsFiles(coreCoords).name;
            splitName = strsplit(coreCoordsName, '_');
            iCore = splitName{1};
            fprintf('Sample: %s - core %s Started', sampleName, iCore);

            croppingdata = load( [ cropCoordsFiles(coreCoords).folder filesep coreCoordsName ] );
            rect = croppingdata.rect;
            boundingBox = [rect(1:2), rect(3:4) - rect(1:2)];
            core = zeros(boundingBox(4), boundingBox(3), numChannelNames);
            box = num2cell(uint16(boundingBox));

            for iChan=1:length(chanelSelected)
                chan=chanelSelected(iChan);
                % Crop core from ome.tif
                core(:,:,iChan) = bfGetPlane(sampleImage, chan, box{:});
            end
            
            %cc = uint32(core);
            cc = uint16(core);
            try
            bfsave(cc,char(strcat(coresFolder, filesep, 'core', iCore, '.tif')), 'dimensionOrder', 'XYCZT');
            catch
                 fprintf('Error detected');
            end     
        end
        toc
end