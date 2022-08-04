function data_struct=HR_fix2(data_struct)
fs=1000;
data_struct=HR_fix(data_struct);
ox_data=double(data_struct.pulse_ox_raw);

ox_data2=ox_data-movmean(ox_data,200);
ox_data2=ox_data2./movmin((movmax(abs(ox_data2),500)),500);
[buffer2,hilb]=BandFilt_Order(ox_data2,fs,fs/2,5,15);
resid=(movmin(movmax(sqrt(movmean((ox_data2-buffer2).^2,500)),1000),1000));
resid=resid-median(resid);
errors=logical(conv(resid>.03,ones(1,500),'same'));

[peaks, troughs,p_amp,t_amp]=findpeaks(buffer2(51:end));
% p_amp_full=interp1(peaks,p_amp,1:length(ox_data),'linear','extrap');
% p_amp_full(errors)=interp1(find(~errors),medfilt1(p_amp_full(~errors),1000),find(errors),'linear','extrap');
% pa_mm=movmean(p_amp_full,1000);
% lower_bound=pa_mm-3*nanstd(p_amp_full./pa_mm).*pa_mm;
% binary_p=nan(1,length(ox_data));
% binary_p(peaks)=p_amp;
% binary_p(binary_p<lower_bound)=NaN;
% peaks=find(~isnan(binary_p));
% p_amp=binary_p(peaks);

% peaks=peaks(p_amp>=5);
peaks=peaks+50;
peaks(find(diff(peaks)<10)+1)=[];
% peaks(ismember(peaks,find(errors)))=NaN;
HR=1e6./[double(data_struct.timestamps(peaks(5:end)))-double(data_struct.timestamps(peaks(1:end-4)))]*60*4;
HR(HR<300 | HR>900)=NaN;
peaktimes=peaks(3:end-2)';
binary_mat=zeros(1,length(data_struct.timestamps));
tofix=ismember(peaktimes,find(errors));
binary_mat(peaktimes(~tofix))=1;
toremove=conv(binary_mat,ones(1,2000),'same')<12;
binary_mat(toremove)=0;
% tofix=conv(tofix,ones(7,1)/7,'same')>0;
% HR(tofix)=NaN;
% tofix=isnan(HR);
peaktimes_good=find(binary_mat>0);
tofix=~ismember(peaktimes,peaktimes_good) | isnan(HR);
HR_full=interp1(peaktimes(~tofix),HR(~tofix),1:length(ox_data),'linear');

peaktimes_good=peaktimes(~tofix);

binary_mat=nan(1,length(data_struct.timestamps));
% peaktimes(tofix)=NaN;
mid_point_pt=min(max(round((peaktimes(2:end)+peaktimes(1:end-1))/2),1),length(binary_mat));
binary_mat(mid_point_pt)=diff(double(data_struct.timestamps(peaktimes))/fs);
tofix2=abs(HR_full-medfilt1(HR_full,3000,'omitnan'))>45;

HR_full(tofix2(:) | errors(:) | toremove(:))=NaN;

HR_full(HR_full<300)=NaN;
binary_mat(isnan(HR_full))=NaN;
pt_new=find(~isnan(binary_mat));
recalc=binary_mat(pt_new);
mid_point_pt=min(max(round((pt_new(2:end)+pt_new(1:end-1))/2),1),length(binary_mat));
binary_mat2=nan(1,length(data_struct.timestamps));
binary_mat2(mid_point_pt)=diff(recalc);

% peaktimes(ismember(peaktimes,find(tofix)))=NaN;
% HR_full=interp1(find(~tofix),HR_full(~tofix),1:length(ox_data),'linear');
HR_full2=HR_full;
HR_full_interp=medfilt1(HR_full,500,'omitnan');
HR_full2(isnan(HR_full))=interp1(find(~isnan(HR_full)),HR_full_interp(~isnan(HR_full)),find(isnan(HR_full)),'linear');
HR_full_interp=LowFilt_Order(HR_full2,1000,10000,.5);
HR_full2(isnan(HR_full))=HR_full_interp(isnan(HR_full));
HR_full2(isnan(HR_full2))=interp1(find(~isnan(HR_full2)),HR_full2(~isnan(HR_full2)),find(isnan(HR_full2)),'nearest','extrap');
peaks_mat=binary_mat2;
% peaks_mat=peaks_mat./nanmean(peaks_mat);
HRV=sqrt(movmean(peaks_mat.^2,1000,'omitnan'));

HRV_interp=medfilt1(HRV,3000,'omitnan');
HRV(isnan(HRV))=interp1(find(~isnan(HRV)),HRV_interp(~isnan(HRV)),find(isnan(HRV)),'linear');
HRV(isnan(HRV))=interp1(find(~isnan(HRV)),HRV(~isnan(HRV)),find(isnan(HRV)),'nearest','extrap');
HRV=movmean(HRV,300);
data_struct.HR_fix=HR_full;
data_struct.HR_fix2=HR_full2;
data_struct.HRV=HRV;
data_struct.HRpeak_vector=peaks_mat;