function [peaks, troughs, peaks_amp, troughs_amp]=findpeaks(x)
x=x(:)';
nantimes=find(isnan(x));
differential=[diff([0 diff(x)]>0) 0];
differential(1)=0;
differential(nantimes)=0;
differential(max(nantimes-1,1))=0;
differential(min(nantimes+1,length(x)))=0;
peaks=find(differential==-1);
troughs=find(differential==1);
len=min(length(peaks),length(troughs));
if len>0
if peaks(1)>troughs(1)
peaks_amp=x(peaks(1:len))-x(troughs(1:len));
troughs_amp=[NaN x(troughs(2:len))-x(peaks(1:len-1))];
else
peaks_amp=[NaN x(peaks(2:len))-x(troughs(1:len-1))];
troughs_amp=x(troughs(1:len))-x(peaks(1:len));
end
if length(troughs_amp)<length(troughs)
    troughs_amp=[troughs_amp NaN];
end
if length(peaks_amp)<length(peaks)
    peaks_amp=[peaks_amp NaN];
end
else
    peaks_amp=[];
    troughs_amp=[];
end