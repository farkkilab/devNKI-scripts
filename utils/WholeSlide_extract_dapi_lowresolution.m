%Script to get the first plane and all the DAPI channels in the folder
%/ fernpere - fperez

basePath = 'D:\users\fperez\NKI_TMAs_AF\';
outputsubfolder = 'slides_dapi'; %Here are the whole slides
omePath = 'registration';
omeSuffix = '.ome.tif';

%Select samples
sampleList = dir( [ basePath 'TMA*' ] );

%Background and DNA channels selection
%DNAchannels=[ 1 5 2 4 ];

DNAchannels=1;

for sample = 1:length(sampleList)
        sampleName = sampleList(sample).name;
        disp(sampleName)
        tic
        %Read ometif image for that sample
        sampleImage =  bfGetReader( [ basePath sampleName filesep omePath filesep sampleName omeSuffix ]);
        %sampleImage.setSeries(4)
        numChannels = sampleImage.getImageCount();
        
        outputFolder = [basePath filesep sampleName filesep outputsubfolder ];
        mkdir(outputFolder);

        for iChan=1:length(DNAchannels)
                chan=DNAchannels(iChan);
                iPlane = sampleImage.getIndex(0, chan - 1, 0) + 1;
                fprintf('Proccesing sample: %s - channel %d \n', sampleName, chan);
                try
                   slide=bfGetPlane(sampleImage, iPlane);
                catch
                    fprintf('Error detected')
                end
                cc = uint16(slide);
                bfsave(cc, char(strcat(outputFolder, filesep, sampleName,'_chanel', int2str(chan), '.tif')));

        end

end