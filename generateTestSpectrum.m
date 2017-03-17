% Generate test spectrum
x=1:1715;

testSpectrum=zeros(1,1715);

% 20 random numbers from 1 to 1715
peakLocations=randi(1715,1,20);

% Random peak amplitudes
peakAmplitudes=randi(30,5,20);

for i=1:20
	testSpectrum=testSpectrum+100/sqrt(peakAmplitudes(i))*exp(-(x-peakLocations(i)).^2/peakAmplitudes(i)^2);
end

plot(x,testSpectrum)