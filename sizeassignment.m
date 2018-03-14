%loading files
load('momentumMatrix.mat');
load('booktomarket.mat');
load('SIZE.mat');
load('monthlyReturns.mat');
%formatting NaN 
Size(Size==-99.9900)=NaN;
booktomarket(booktomarket==-99.9900)=NaN;
momentumMatrix(momentumMatrix==-100)=NaN;
monthlyReturns1(monthlyReturns1==-99.9900)=NaN;
%formating table sizes
Size=Size(1:1098,:);
Size=Size(13:1098,:);
monthlyReturns1 = monthlyReturns1(13:1098,:);
%summary statistics table 1
n = size(Size,1);
cross = size(Size,2);
CrossSectionalMean = nanmean(Size,2);
MeanThroughTime = mean(CrossSectionalMean);
StandardDeviation = std(CrossSectionalMean);
Skewness = skewness(CrossSectionalMean);
Kurtosis = kurtosis(CrossSectionalMean);
percentile5 = prctile(CrossSectionalMean,5);
percentile25 = prctile(CrossSectionalMean,25);
Median = median(CrossSectionalMean);
percentile75 = prctile(CrossSectionalMean,75);
percentile95 = prctile(CrossSectionalMean,95);
maximum = max(CrossSectionalMean);
minimum = min(CrossSectionalMean);
%end table 1

%table 2 correlations 
%momentum cross sectional 
Mom = nanmean(momentumMatrix,2); 
%Pearson Correlation 
pearson = corr(CrossSectionalMean,Mom);
%Sprearman Correlation
spearman = corr(CrossSectionalMean,Mom,'type', 'Spearman');
%end table 2

%table 3 panel A
treatNan = isnan(Size)|isnan(monthlyReturns1);
sizeMat = zeros(10,size(Size,1));
for i=1:size(Size,1)
    [sorted,index] = sort(Size(i,~treatNan(i,:)),2);
    NumofAssetMonth = size(sorted,2);
    sizeMat(1,i) = mean(sorted(1:round(0.1*NumofAssetMonth)));
    sizeMat(2,i) = mean(sorted(1:round(0.1*NumofAssetMonth)+1 : round(0.2*NumofAssetMonth)));
    sizeMat(3,i) = mean(sorted(1:round(0.2*NumofAssetMonth)+1 : round(0.3*NumofAssetMonth)));
    sizeMat(4,i) = mean(sorted(1:round(0.3*NumofAssetMonth)+1 : round(0.4*NumofAssetMonth)));
    sizeMat(5,i) = mean(sorted(1:round(0.4*NumofAssetMonth)+1 : round(0.5*NumofAssetMonth)));
    sizeMat(6,i) = mean(sorted(1:round(0.5*NumofAssetMonth)+1 : round(0.6*NumofAssetMonth)));
    sizeMat(7,i) = mean(sorted(1:round(0.6*NumofAssetMonth)+1 : round(0.7*NumofAssetMonth)));
    sizeMat(8,i) = mean(sorted(1:round(0.7*NumofAssetMonth)+1 : round(0.8*NumofAssetMonth)));
    sizeMat(9,i) = mean(sorted(1:round(0.8*NumofAssetMonth)+1 : round(0.9*NumofAssetMonth)));
    sizeMat(10,i) = mean(sorted(1:round(0.9*NumofAssetMonth)+1 :end));
end 
average=mean(sizeMat,2);
  
% Table 3 panel b 
%soritng montlhy returns according to size 
[sorted,index] = sort(Size,2);
for i=1:size(Size,1)  
   for j=1:size(Size,2)
       Returns(i,j) = monthlyReturns1(i,index(i,j));
   end
end

%now the Nan values are in the last rows in each line in matrix Returns
treatNan = isnan(Returns);
%estimate the mean return in each decile every month
EquallyMeanDecileMatrix = zeros(1086,10);
WeightMeanDecileMatrix = zeros(1086,10);
for i=1:size(Returns,1) 
    NumofAssetMonth =  size(Returns,2) - sum(treatNan(i,:));
    NumofAssetInDecileRound = round( 0.1 * NumofAssetMonth);
    EquallyMeanDecileMatrix(i,1) = mean(Returns(i,1:round(0.1*NumofAssetMonth)));    
    t=0.1;
    l=0.2;
    for j=2:9
        EquallyMeanDecileMatrix(i,j) = mean(Returns(i,round(t*NumofAssetMonth)+1 : round(l*NumofAssetMonth)));
        t=t+0.1; 
        l=l+0.1;    
    end
    EquallyMeanDecileMatrix(i,10) = mean(Returns(i,round(0.9*NumofAssetMonth)+1 : NumofAssetMonth));
end

%table 2 Panel A
for i=1:10
MeanDecilethroughTime(1,i) = mean(EquallyMeanDecileMatrix(:,i));
end

WeightMeanDecileMatrix = zeros(1086,10);
WeightMatrix = Returns .* sorted ;    
for i=1:size(Returns,1) 
    NumofAssetMonth =  size(Returns,2) - sum(treatNan(i,:));
    NumofAssetInDecileRound = round( 0.1 * NumofAssetMonth);
    WeightMeanDecileMatrix(i,1) = sum(WeightMatrix(i,1:round(0.1*NumofAssetMonth))) / sum(sorted(i,1:round(0.1*NumofAssetMonth)));   
    t=0.1;
    l=0.2;
    for j=2:9
        WeightMeanDecileMatrix(i,j) = sum(WeightMatrix(i,round(t*NumofAssetMonth)+1 : round(l*NumofAssetMonth))) / sum(sorted(i,round(t*NumofAssetMonth)+1 : round(l*NumofAssetMonth)));
        t=t+0.1; 
        l=l+0.1;    
    end
    WeightMeanDecileMatrix(i,10) = sum(WeightMatrix(i,round(0.9*NumofAssetMonth)+1 : NumofAssetMonth)) / sum(sorted(i,round(0.9*NumofAssetMonth)+1 : NumofAssetMonth));
end
for i=1:10
WeightMeanDecilethroughTime(1,i) = mean(WeightMeanDecileMatrix(:,i));
end
% Montlhy Returns of equally balanced and weight balanced Size portfolio 
SizeMonthlyReturn = EquallyMeanDecileMatrix(:,1) - EquallyMeanDecileMatrix(:,10);
WeightSizeMonthlyReturns =  WeightMeanDecileMatrix(:,1) -  WeightMeanDecileMatrix(:,10);   
save SizeMonthlyReturn
save WeightSizeMonthlyReturns

%Regressions 
%loading of ff factors 
load('FFResearchDataFactors.mat'); 
FFFactors = FFResearchDataFactors(13:1098,:);

MarketPremium = FFFactors(:,1) - FFFactors(:,4);
SMB = FFFactors(:,2);
HML = FFFactors(:,3);
Y = [MarketPremium, SMB, HML];

EquallyWeightedResultsCAPM = regstats(SizeMonthlyReturn,MarketPremium);
EquallyWeightedResultsFFF = regstats(SizeMonthlyReturn,Y);

SizeWeightedResultsCAPM = regstats(WeightSizeMonthlyReturns,MarketPremium);
SizeWeightedResultsFFF = regstats(WeightSizeMonthlyReturns,Y);

%compute nwse for eq
EqualCoefCAPM = EquallyWeightedResultsCAPM.tstat.beta;
EqualResidaulCAPM = EquallyWeightedResultsCAPM.r;
newstandarerrorsEqCAPM = nwse(EqualResidaulCAPM,MarketPremium);
tstatCAPMEq = EqualCoefCAPM ./ newstandarerrorsEqCAPM;
save tstatCAPMEq;

EqualCoefFFF = EquallyWeightedResultsFFF .tstat.beta;
EqualResidaulFFF = EquallyWeightedResultsFFF .r;
newstandarerrorsEq = nwse(EqualResidaulFFF,Y);
tstatFFFEq = EqualCoefFFF ./ newstandarerrorsEq;
save tstatFFFEq;

%compute nwse for size weighted
WeightCoefCAPM = SizeWeightedResultsCAPM.tstat.beta;
WeightResidaulCAPM = SizeWeightedResultsCAPM.r;
newstandarerrorsWeightCAPM = nwse(WeightResidaulCAPM,MarketPremium);
tstatCAPMWeight = WeightCoefCAPM ./ newstandarerrorsWeightCAPM;
save tstatCAPMWeight;

WeightCoefFFF = SizeWeightedResultsFFF.tstat.beta;
WeightResidaulFFF = SizeWeightedResultsFFF.r;
newstandarerrorsWeightFFF = nwse(WeightResidaulFFF,Y);
tstatFFFWeight = WeightCoefFFF ./ newstandarerrorsWeightFFF;
save tstatFFFWeight;

%end 
