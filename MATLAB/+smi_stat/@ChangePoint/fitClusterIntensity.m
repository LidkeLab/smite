function icp=fitClusterIntensity( data, logBayesThreshold, plotit)
    if nargin==2
        plotit=false;
    end
    data=round(data);
    firstZero=find(~data,1);
    if ~isempty(firstZero)
        assert(firstZero>1); %The data does not start with a 0
        positiveData=data(1:firstZero-1);
        assert(all(positiveData>0)); %Check we got all non-zero values
        assert(all(data(firstZero:end)==0)); %Check we did not leave any photons behind
        data=positiveData;
    end
    data=data(3:end);
    icp=IntensityCPA(data, logBayesThreshold);
    if plotit
        icp.plotIntensityEstimate();
    end
end
