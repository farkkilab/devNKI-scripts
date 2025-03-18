% Crop TMA ome.tif into the already detected cores
% Quantify with already computed masks.
% /casado

basePath = 'D:\users\fperez\NKI_TMAs_AF\';
maskPath = 'whole-cell-segmentation2'; %Input folder name
maskFileName = '_Probabilities_cytomask2.tiff';
omePath = 'registration';
omeSuffix = '.ome.tif';
outputsubfolder = 'quantification2';
cropCoordsPath = 'dearray\cropCoords\';
cropCoordsFileName = '*_cropCoords.mat';

channelNames = readtable( [basePath filesep 'channel_list.csv'], 'ReadVariableNames', false);
numChannelNames = size(channelNames, 1); % this is only to allocate memory before loops

% Quantification features
%  - The eccentricity is the ratio of the distance between the foci of 
%    the ellipse and its major axis length. (0 = circle, 1 = segment).
%  - Solidity: Proportion of the pixels in the convex hull that are also in the region
%    . Computed as Area/ConvexArea. Solidity is useful to quantify the 
%    amount and size of concavities in an object boundary. Holes are also 
%    often included. For example, it distinguishes a star from a circle,
%    but doesn?t distinguish a triangle from a circle.
%  - Roundness: C = 4 * pi Area / Perimeter^2 (to replace Solidity in the
%    next datasets.
featureNames = {'Area', 'Eccentricity', 'Perimeter', 'Solidity', 'MajorAxisLength', 'MinorAxisLength'};
XYnames = {'X_position','Y_position'};

% List of samples
sampleList = dir( [ basePath 'TMA*' ] );

selected = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
%selected = 10;

parfor sample = 1:length(selected)
%for sample = 1
    samp = selected(sample);
    sampleName = sampleList(samp).name;
    disp(sampleName)
    tic
    
    %Read large ome.tif
    sampleImage =  bfGetReader( [ basePath sampleName filesep omePath filesep sampleName omeSuffix ] );
    numChannels = sampleImage.getImageCount();
    
    if numChannelNames ~= numChannels
        disp('ERROR: number of channel names and actual channels do not match');
        continue
    end
    % List available coordinate files
    cropCoordsFiles = dir( [ basePath filesep sampleName filesep cropCoordsPath filesep cropCoordsFileName ] );
    
    % Create a folder for output core ome.tif files in dearray folder
    %coresFolder = [ basePath filesep sampleName filesep cropCoordsPath filesep 'cores' ];
    %mkdir(coresFolder);
    
    % Create a folder for quantification step
    outputFolder = [basePath filesep sampleName filesep outputsubfolder ];
    % init alldata table to store all cores together, 'cellid channels
    % morphFeatures roundness x y'
    alldata = [];
    mkdir([ outputFolder filesep 'cores'] );
    
    for coreCoords = 1:length(cropCoordsFiles)
    %for coreCoords=2
        coreCoordsName = cropCoordsFiles(coreCoords).name;
        splitName = strsplit(coreCoordsName, '_');
        iCore = splitName{1};
        
        fprintf('Sample: %s - core %s Started \n', sampleName, iCore);
        % Coordinate .mat files must contain a 'rect' object
        croppingdata = load( [ cropCoordsFiles(coreCoords).folder filesep coreCoordsName ] );
        rect = croppingdata.rect;
        % Load mask
        mask = imread( [ basePath sampleName filesep maskPath filesep 'core' iCore maskFileName ] );
        l = size(mask);
        if (l(1) == rect(4) - rect(2))
            boundingBox = [rect(1:2), rect(3:4) - rect(1:2)];
        else
            boundingBox = [rect(1:2), rect(3:4) - rect(1:2) + 1 ];
        end
               
        core = zeros(boundingBox(4), boundingBox(3), numChannels); 
        if(boundingBox(3) ~= boundingBox(4))
            disp('First edge core, check it out')
        end
        box = num2cell(uint16(boundingBox));
        % Parfor doesn't allow initializing the bioformats reader out of
        % the loop
        sampleImage =  bfGetReader( [ basePath sampleName filesep omePath filesep sampleName omeSuffix ] );
         % Iterate over channels
        for iChan=1:numChannelNames
            % Crop core from ome.tif
            core(:,:,iChan) = bfGetPlane(sampleImage, iChan, box{:});
        end
        
        % save core ome.tif
        %bfsave(core,char(strcat(coresFolder, filesep, 'core', iCore, '.tif')),'BigTiff',true);
        
        % Quantify channels
        getMeanFunction = @(iChan) struct2array(regionprops(mask, core(:,:,iChan), 'MeanIntensity'))';
        meanIntensities = cell2mat(arrayfun(getMeanFunction,1:numChannelNames, 'UniformOutput',0));

        % Calculate morphological features
        shapeFeatures = regionprops(mask, featureNames);
        roundness = (4 * pi * [shapeFeatures.Area]') ./ [shapeFeatures.Perimeter]'.^2;
        % Calculate X and Y positions
        xystruct = regionprops(mask, 'Centroid');
        xytemp = cat(1, xystruct.Centroid);
        
        % Write all variables as CSV file
        lenIDs = unique(mask);
        CellId = array2table(double(lenIDs(lenIDs ~= 0)), 'VariableNames', {'CellId'});
        CoreId = array2table( repmat(string(iCore), size(CellId,1), 1), 'VariableNames', {'CoreId'});
        meanData = array2table(meanIntensities,'VariableNames', table2cell(channelNames) );
        morphData = struct2table(shapeFeatures);
        roundData = array2table(roundness, 'VariableNames',{'Roundness'} );
        xyData = array2table( [xytemp(:,1) xytemp(:,2) ], 'VariableNames', XYnames);
        
        %Next is to detect jumps of pixels values in the mask 
        %A jump is when a instensity value is skiped from the mask (1,2,4,5), so that
        %cell (3) dont exists
        %When there is a jum, that row will contain only NaN+
        [rows, columns] = size(meanData);
        for row=1:rows
            if (all(isnan(meanIntensities(row,:))))
                meanData(row,:) = [];
                morphData(row,:) = [];
                roundData(row,:) = [];
                xyData(row,:) = [];
            end
        end
        
        output = [ CoreId, CellId, meanData, morphData, roundData, xyData ];
        
        writetable( output, [ outputFolder filesep 'cores' filesep iCore] , 'Delimiter', '\t');
        alldata = [alldata ; output ];
    end
    writetable(alldata, [outputFolder filesep sampleName '.csv'], 'Delimiter', '\t');
    toc
end