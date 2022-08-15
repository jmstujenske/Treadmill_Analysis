function [Filtered,FiltHilb,FiltAmp,FiltPh] = HighFilt_Order(eeg,sampFreq,order,high)

Nyquist = sampFreq/2;
MyFilt=fir1(order,[high]/Nyquist,'high');

Filtered = Filter0(MyFilt,eeg);

if (nargout>1)
    FiltHilb = hilbert(Filtered);
    FiltPh = angle(FiltHilb);
    FiltAmp = abs(FiltHilb);
end
