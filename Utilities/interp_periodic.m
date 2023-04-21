function Y=interp_periodic(X,fr,filt_cut,offset)
%interp_periodic(X,fr,filt_cut)
%Intelligent interpolation of known periodic signal
%
%IN:
%X - signal (vector) with NaN for values to interpolate
%fr - sampling rate
%filt_cut - filter cutoffs
%   [0 f] = low-pass
%   [f inf] = high-pass
%   [f1 f2] = band-pass
%offset - number of angles on either side to ignore for interpolation
%
%OUT:
%Y - signal with interpolated values replaced for NaNs
%
%This method breaks the signal into a periodic and non-periodic component
%using the Hilbert transform. The periodic component is interpolated
%separately from the non-periodic component (which itself is often related
%to the periodic component so its own periodicity is taken ito account).
%The components are then recombined to get the final interpolation.
%In a last step, due to rare discontinuities of the interpolation with the
%underlying signal, discontinuities are fixed by rescaling individual
%segments
%
if nargin<4 || isempty(offset)
    offset=0;
end
if nargin<3
    filt_cut=[];
end
%initially interpolate with a spline; this will do a poor job, but we need to
%remove the NaNs for the nextsteps
tofix=isnan(X); %mark the original missing points as "tofix"
X(tofix)=interp1(find(~tofix),X(~tofix),find(tofix),'spline','extrap');

Y=X; %Y is a copy of the original input, in case we want to make any changes to X
% X=sgolayfilt(X,1,25); %optional smoothing

%
if isempty(filt_cut)
    hilb=Hilbert(X);
    filt=X;
else
    if length(filt_cut)<2
        warning('One input for filt_cut is ambiguous, will do a highpass.');
        filt_cut=[filt_cut inf];
    end
    [filt,hilb]=customfilter(X,fr,filt_cut);
end
angle_interpolated=interpolate_angles(hilb,tofix,offset);

%have to scale the periodic signal to make it match up with actual trace,
%and will need to interpolate this scaling
%The power envelope is also periodic and related to the period of the
%underlying signal
power_trace=(filt./cos(angle_interpolated));

tofix2=logical(conv(tofix,ones(1,ceil(fr*.5/2)),'same')); %extend periods needing fixing by 250ms
power_trace(tofix2)=NaN;

power_trace_interp=power_trace;
power_trace_interp(tofix2)=interp1(find(~tofix2),power_trace(~tofix2),find(tofix2),'linear','extrap');

[filt2,hilb2]=customfilter(power_trace_interp,fr,filt_cut);
angle_interpolated2=interpolate_angles(hilb2,tofix,offset);

%determine relationship between the power envelope and periodicity of the
%signal by regression
angle_angle_rel=MeanSmooth_3D_circular([angle_interpolated],[angle_interpolated2],-pi-pi/40:pi/20:pi+pi/40);

angle_interpolated2_fixed=interp1(-pi:pi/20:pi,angle_angle_rel(1:end-1),cos(angle_interpolated));
angle_interpolated2_fixed=(angle_interpolated2_fixed-nanmin(angle_interpolated2_fixed))/2;

baseline=movmax(movmin(power_trace,ceil(fr),'omitnan'),ceil(fr),'omitnan');
power_trace_fix=angle_interpolated2_fixed+baseline;
power_trace_fix=(power_trace_fix-baseline).*movmin(movmax(power_trace-baseline,fr,'omitnan'),fr)+baseline;
power_trace(tofix2)=power_trace_fix(tofix2);

%residual between original signal and the filtered version
newfilt=power_trace_fix.*cos(angle_interpolated);
resid=X-newfilt;
resid(tofix)=interp1(find(~tofix),resid(~tofix),find(tofix),'linear');

%finally, interpolate the signal using the rescaled periodic signal
Y(tofix)=power_trace(tofix).*cos(angle_interpolated(tofix))+resid(tofix);
Y_all=newfilt+resid;

%Lastly, sometimes the interpolation is slightly off at the ends, so will
%linearly rescale each interpolation to match 
ins=find(tofix);
ends=find(diff(ins)>1);
starts=[ins(1);min(ins(ends+1),length(tofix))];
ends=[ins(ends);ins(end)];

for a=1:length(starts)
    initdiff=Y_all(max(starts(a)-1,1))-Y(max(starts(a)-1,1));
    lastdiff=Y_all(min(ends(a)+1,length(Y_all)))-Y(min(ends(a)+1,length(Y_all)));
    slope=linspace(initdiff,lastdiff,(ends(a)+3-starts(a)))';
    Y(starts(a):ends(a))=Y(starts(a):ends(a))-slope(2:end-1);
end

function angle_interpolated=interpolate_angles(hilb,tofix,offset)
%angle_interpolated will be the angles of the hilbert while angle_trace will
%be the "predicted angle" based on frame-by-frame changes in the hilbert.
%These will end up out of phase due to mistakes in the prediction, but
%they should be the same frequency

angle_interpolated=angle(hilb);
angle_trace=angle_interpolated;
qs=quantile(mod(angle_trace+pi,2*pi),0:.01:1);
angle_sym=interp1(qs,0:.01:1,mod(angle_trace+pi,2*pi))*2*pi-pi;
breaks=find(diff(angle_sym)>pi);
for a=1:length(breaks)
    angle_sym(breaks(a)+1:end)=angle_sym(breaks(a)+1:end)-2*pi;
end
breaks=find(diff(angle_sym)<-pi);
for a=1:length(breaks)
    angle_sym(breaks(a)+1:end)=angle_sym(breaks(a)+1:end)+2*pi;
end
angle_sym=[0;diff(angle_sym)];
angle_sym(tofix)=NaN;

%get rid of spurious angle changes
angle_sym=movmedian(angle_sym,10,'omitnan');

%need to locally average because the respiration trace is asymmetric, and
%therefore the frame-by-frame angle trace is actually periodic itself
angle_sym=movmean(angle_sym,120,'omitnan');

angle_trace=angle_sym;

%loop through missing segments
%in each loop, predict the angle change vs the actual angle change
%if off by more than half of a cycle [modifiable], then add full cycles to
%make up for the difference
angle_trace(isnan(angle_trace))=0;
angle_trace=cumsum(angle_trace);
ins=find(tofix);
ends=find(diff(ins)>1);
starts=[ins(1);min(ins(ends+1),length(tofix))];
ends=[ins(ends);ins(end)];

for a=1:length(starts)
    expectedangletraverse=diff(angle_trace([max(starts(a)-1-offset,1) min(ends(a)+1+offset,length(angle_trace))]));
    angle_start=angle(hilb(max(starts(a)-1-offset,1)));
    angle_end=angle(hilb(min(ends(a)+1+offset,length(angle_trace))));
    angle_traverse=circ_dist(angle_end,angle_start);
    if angle_traverse<0
        angle_traverse=angle_traverse+2*pi;
    end
    extracycles=(expectedangletraverse-angle_traverse)/(2*pi);
    if mod(extracycles,1)>.5
        extracycles=ceil(extracycles);
    else
        extracycles=floor(extracycles);
    end
    extracycles=max(extracycles,0);
    totalangle=extracycles*2*pi+angle_traverse;
    if starts(a)<2+offset
        modifier=1+offset-starts(a)+1;
    else
        modifier=0;
    end
    newangles=linspace(angle_start+pi,angle_start+totalangle+pi,ends(a)-starts(a)+3-modifier+2*offset);
    newangles=mod(newangles,2*pi)-pi;
    reverse_sym=interp1((0:.01:1)*2*pi-pi,qs-pi,newangles,'linear','extrap');
    angle_fix=linspace(reverse_sym(1)-newangles(1),reverse_sym(end)-newangles(end),ends(a)-starts(a)+3-modifier+2*offset);
    angle_interpolated(max(starts(a)-1-offset,1):ends(a)+1+offset)=reverse_sym-angle_fix;
end
angle_interpolated=mod(angle_interpolated+pi,2*pi)-pi;
angle_interpolated=angle_interpolated(1:length(hilb));

function [filt,hilb]=customfilter(X,fr,filt_cut)
if ~isempty(filt_cut)
    if isinf(filt_cut(2))
        [filt,hilb]=HighFilt_Order(X,fr,fr/2,filt_cut(1));
    elseif filt_cut(1)==0
        [filt,hilb]=LowFilt_Order(X,fr,fr/2,filt_cut(2));
    elseif ~isempty(filt_cut)
        [filt,hilb]=BandFilt_Order(X,fr,fr/2,filt_cut(1),filt_cut(2));   
    end
else
    filt=[];
    hilb=[];
end
