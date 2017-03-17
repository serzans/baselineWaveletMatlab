% Load example spectrum or test spectrum
STest=transpose(importdata('ExampleSpectra.txt'));
save exampleScpetrum.mat STest
% STest=testSpectrum;


% Array containing the signal S (initialise as GPU array for faster
% calculations)
S=gpuArray(STest);

% Maximum scale factor for the CWT
maxA=70;

% Generate Mexican hat wavelet CWT coefficient matrix
mexicanCStruct=cwtft(S,'scales',{1,1,maxA,'lin'},'wavelet','mexh');
mexicanC=real(mexicanCStruct.cfs);

% Save array to binary file
save mexicanC.mat mexicanC

% Find local maxima in every row
maxC=findLocalMaxima(mexicanC);

% Find connected regions in the binary image
% Accumulate a list of ridges with pixel coordinates
ridgesUnfiltered=regionprops(bwconncomp(maxC,8),'pixellist'); % Gives fragmented ridges

% For plotting the ridges
ridgeBinaryMatrix=zeros(size(mexicanC));


% Get statistics on connected regions (ridges)
for j=1:size(ridgesUnfiltered)

	% For plotting the ridges...
	currentRidge=fliplr(ridgesUnfiltered(j).PixelList);
	for i=1:size(currentRidge)
		x=currentRidge(i,1);
		y=currentRidge(i,2);
		ridgeBinaryMatrix(x,y)=1;
	end

	% Get statistics on ridges. Some unused parameters commented out.
	[params]=findRidgeParameters(ridgesUnfiltered(j).PixelList,mexicanC);
% 	ridgesUnfiltered(j).verticalLength=params(1);
% 	ridgesUnfiltered(j).horizontalDeviation=params(2);
	ridgesUnfiltered(j).pixelCount=params(3);
% 	ridgesUnfiltered(j).coefficientSum=params(4);
% 	ridgesUnfiltered(j).interceptChannel=params(5);
    ridgesUnfiltered(j).maxA=params(6);
%     ridgesUnfiltered(j).maxCoordinate=params(7);
    ridgesUnfiltered(j).maxValue=params(9);
end


% Accumulate filtered ridges in the cell struct thresholdedRidges
thresholdedRidges={};

% Output matrix for plotting filtered ridges (important for dilation)
ridgeBinaryMatrixFiltered=zeros(size(mexicanC));

for i=1:size(ridgesUnfiltered)
   % Culling ridges containing few pixels near small scall factors 'a'
   minimumPixelCount=3;
   maxScaleFactorThreshold=10;
   minCoefficientThreshold=2500;
   
   if ridgesUnfiltered(i).maxA<maxScaleFactorThreshold
       if ridgesUnfiltered(i).pixelCount<minimumPixelCount 
       end
   % Selecting for max coefficient value on the ridge
   elseif ridgesUnfiltered(i).maxValue<minCoefficientThreshold
   % Pass
   else
        thresholdedRidges{end+1}=ridgesUnfiltered(i);
        
        % For plotting
       	currentRidge=fliplr(ridgesUnfiltered(i).PixelList);
        for k=1:size(currentRidge)
            x=currentRidge(k,1);
            y=currentRidge(k,2);
            ridgeBinaryMatrixFiltered(x,y)=1;
        end
   end 

end

% Increase ridge thickness to connect fragmented regions
ridgeBinaryMatrixFiltered=imdilate(ridgeBinaryMatrixFiltered,strel('disk',4));
% Decrease ridge thickness
ridgeBinaryMatrixFiltered=imerode(ridgeBinaryMatrixFiltered,strel('disk',2));

% Find connected regions in the filtered binary image
% Accumulate a list of ridges with pixel coordinates
ridgesFiltered=regionprops(bwconncomp(ridgeBinaryMatrixFiltered,8),'pixellist');



% Get statistics on the filtered connected regions
for j=1:size(ridgesFiltered,1)
	params=findRidgeParameters(ridgesFiltered(j).PixelList,mexicanC);
	ridgesFiltered(j).verticalLength=params(1);
	ridgesFiltered(j).horizontalDeviation=params(2);
	ridgesFiltered(j).pixelCount=params(3);
	ridgesFiltered(j).coefficientSum=params(4);
	ridgesFiltered(j).interceptChannel=params(5);
    ridgesFiltered(j).maxA=params(6);
    ridgesFiltered(j).maxChannel=params(7);
    ridgesFiltered(j).maxScaleFactor=params(8);
end


% Peak width estimation...
% Generate Haar wavelet transform coefficient matrix (cannot use gpuArray for this)
haarAbsoluteC=abs(cwt(gather(S),1:maxA,'haar'));

% Save data to binary file
save haarAdbsoluteC.mat haarAbsoluteC


% Find all local minima
minimaBinaryMatrix=zeros(size(mexicanC));

% Iterate over rows in the Haar wavelet CWT & find minima
for i=1:maxA
    % Use findpeaks function with the negative to find minima
    [pks,locs]=findpeaks(-1*haarAbsoluteC(i,:));
    minimaBinaryMatrix(i,locs)=1;
end


% Combine minima with location of ridges to get a logical mask
combinedMatrix=or(ridgeBinaryMatrixFiltered,minimaBinaryMatrix);


% For each of the peak maxima find the closest mask on the left and on the
% right. Estimate the peak width from the difference.
for j=1:size(ridgesFiltered,1)
    maxChannel=ridgesFiltered(j).maxChannel;
    maxScaleFactor=ridgesFiltered(j).maxScaleFactor;
    tempMatrix=combinedMatrix;
    
    % Mask out the current ridge
    currentRidge=fliplr(ridgesFiltered(j).PixelList);
        for k=1:size(currentRidge)
            x=currentRidge(k,1);
            y=currentRidge(k,2);
            tempMatrix(x,y)=0;
        end
    
    % Find distances in the row to maximum coordinate (sorted)
    distances=sort(abs(find(tempMatrix(maxScaleFactor,:))-maxChannel));
    % Add two closest distances to get peak width
    ridgesFiltered(j).width=distances(1)+distances(2)+1;
end

figure
plot(S)
hold on
for i=1:size(ridgesFiltered,1)
    width=ridgesFiltered(i).width;
    center=ridgesFiltered(i).maxChannel;
    plot([center-width/2,center+width/2],[2,2],'g');
    plot([center,center],[1,11000],'r');
end    


% Plotting
% Consider plotting ridge maxima points
% figure;
% subplot(3,1,1);
% imagesc(mexicanC);
% colorbar
% 
% subplot(3,1,2);
% imagesc(ridgeBinaryMatrix);
% 
% subplot(3,1,3);
% imagesc(ridgeBinaryMatrixFiltered);



%%% Suggested methods for paralelisation!!!!!
% Use gpuarray to generate MHW CWTs
% Use a parfor array to generate the Haar CWTs 

% Use gpuarray for functions bwconncomp, imdilate and imerode
% (finding connected regions and dealing with fragmentation of ridges)

% Use parfor for collecting statistics on ridges (function
% findRidgeParameters)







