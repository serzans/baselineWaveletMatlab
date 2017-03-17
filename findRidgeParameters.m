function [ridgeParameters] = findRidgeParameters(ridgePixelList,C)
	% ridgePixelList - list of pixel coordinates

	% Vertical length
	verticalLength=max(ridgePixelList(:,2))-min(ridgePixelList(:,2))+1;

	% Horizontal deviation
	horizontalDeviation=max(ridgePixelList(:,1))-min(ridgePixelList(:,1))+1;

	% coefficientSum (sum of wavelet coeffs.)
	coefficientValues=zeros(size(ridgePixelList,1),1);
	for i=1:size(ridgePixelList,1)
		coefficientValues(i)=gather(C(ridgePixelList(i,2),ridgePixelList(i,1)));
    end
    coefficientSum=sum(coefficientValues);
    
    % Coordinate of maximum value on the ridge
    [maxValue,idx]=max(coefficientValues);
    maxCoordinate=ridgePixelList(idx,:);

	% Pixel Count
	pixelCount=size(ridgePixelList,1);

	% Channel at intercept a=0
	[MinA,minValueIndex]=min(ridgePixelList(:,2));
    
    % Maximum scale factor 'a' of the ridge
    maxA=max(ridgePixelList(:,2));
    
	% Approximate intercept channel at a=1
	interceptChannel=ridgePixelList(minValueIndex,1);

	ridgeParameters=[verticalLength,horizontalDeviation,pixelCount,coefficientSum,interceptChannel,maxA,maxCoordinate,maxValue];
end