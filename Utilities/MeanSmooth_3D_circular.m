% [m Bins sd stderr n] = MeanSmooth(x,y,Bins)
% splits x into bins, and calculates the mean of y in each bin.
%
% Bins may be a vector in which case the nth bin is all those points greater than
% or equal to x(n) but less than x(n+1), and the last one is all those greater than
% x(end).  Or it may be a number - in which case the range of x is divided into
% bins containing that many points.
%
% third and fourth optional output arguments are the sd and standard error for
% the relevant bin.
%
% see also MedianSmooth

function [m, Bins, sd, stderr, n_count Save] = MeanSmooth_3D_circular(x,y,Bins)

if (length(Bins) == 1)
	PointsPerBin = Bins;
	sorted = sort(x);
	Bins = sorted(1:PointsPerBin:end);
end
nBins = length(Bins);
n_count=zeros(nBins,1);
Save=cell(nBins,1);
m = NaN*ones(nBins,size(y,2));
sd = zeros(nBins,size(y,2));
stderr = zeros(nBins,size(y,2));
for n=1:nBins-1
	PointsInBin = find(x>=Bins(n) & x<Bins(n+1));
	if ~isempty(PointsInBin)
        n_count(n) = length(PointsInBin);
		m(n,:) = circ_mean(y(PointsInBin,:),[],1);
        Save{n}=y(PointsInBin,:);
		sd(n,:) = nanstd(y(PointsInBin,:),1);
		stderr(n,:) = sd(n,:) / sqrt(length(PointsInBin)-1);
	end
end

% do last bin
PointsInBin = find(x>=Bins(nBins));
if ~isempty(PointsInBin)
    Save{nBins}=y(PointsInBin,:);
    n_count(nBins) = length(PointsInBin);
	m(nBins,:) = circ_mean(y(PointsInBin,:),[],1);
	sd(nBins,:) = nanstd(y(PointsInBin,:),1);
	stderr(nBins,:) = sd(nBins,:) / sqrt(length(PointsInBin)-1);
end