function [Filtered,FiltHilb,FiltAmp,FiltPh] = LowFilt_Order(eeg,sampFreq,order,low)

Nyquist = sampFreq/2;
MyFilt=fir1(order,[low]/Nyquist,'low');

Filtered = Filter0(MyFilt,eeg);

if (nargout>1)
    FiltHilb = hilbert(Filtered);
    FiltPh = angle(FiltHilb);
    FiltAmp = abs(FiltHilb);
end
