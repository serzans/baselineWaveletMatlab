function [localMaximaC] = findLocalMaxima(C)
	% Output matrix containing locations of local maxima
	localMaximaC=zeros(size(C));
    
    maxA=size(C,1);
    for a=1:maxA
		currentRow=gather(real(C(a,:)));
        % Find local maxima and their locations
        [pks,pklocs]=findpeaks(currentRow);
        localMaximaC(a,pklocs)=1;
    end
end

