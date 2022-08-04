function data_struct=HR_fix(data_struct)
ox_data=double(data_struct.pulse_ox_raw);
% error_index=movmin(movmax(movmean(abs(ox_data-(max(ox_data)+min(ox_data))/2),200),1000),1000);
% error_index=movmin(movmax(abs((diff(ox_data))),2000),2000);
error_index=(movmean(ox_data>(max(ox_data(1:end-2000))-5) | ox_data<(min(ox_data(1:end-2000))+5),500));
error_index=error_index-movmin(error_index,20000);
ox_data=sgolayfilt(ox_data,3,21);
% ox_data=conv(ox_data,gausskernel(50,10),'same');
% fs=floor((length(data_mat)-1)/diff(data_mat([1 end],1))*1e6);
fs=1000;
% [yo, fo, to]=mtcsg(ox_data./sqrt(movmean(ox_data.^2,20000)),fs*5,fs,fs,fs-500,1.5);
% [yo, fo, to]=mtcsg(ox_data,fs*5,fs,fs,fs-500,1.5);
% ox_data=ox_data-movmean(ox_data,20000);

[buffer2]=BandFilt_Order(ox_data,fs,300,5,15);
% [buffer3]=BandFilt_Order(ox_data,fs,2000,.5,2);

high_amp=movmin(movmax(abs(buffer2),150),150);
% [~,in]=max(yo(fo>2 &fo<=15.5,:).*fo(fo>2 &fo<=15.5));
% to=to*fs+fs/2;

% [buffer]=BandFilt_Order(ox_data,fs,100,30,100);
% low_amp=movmin(movmax(abs(buffer),100),100);
% error_index=sqrt(low_amp);


% error_index=(movmean(nanmean(yo(fo>2 & fo<5,:).*fo(fo>2 & fo<5),1),1));
% error_index=error_index-movmin(error_index,50);
% 
% error_index=sqrt(error_index);
% high_amp=sqrt(high_amp);
% error_index=sqrt(error_index);
% errors=error_index>(nanmean(error_index)+2*nanstd(error_index));
% errors=error_index>340;
% errors=conv(error_index>.03,ones(500,1),'same');
errors=error_index>.03;
signalloss=movmin(movmax((high_amp<10)*1000,1000),1000);
tofix3=conv([false;(double(data_struct.position_linear(3:end))-double(data_struct.position_linear(1:end-2)))>0;false],ones(1000,1),'same')>3;
errors=errors | signalloss |tofix3;
% errors=conv(movmin(movmax(errors*1000,1000),1000),ones(1,500),'same')>0;
errors=(movmin(movmax(errors*1000,1000),1000))>0;

% errors=conv(errors,ones(1,5),'same')>0;
% errors=errors' | fo(in+find(fo<=2,1,'last'))<5;
% yo(:,errors)=[];
% in(errors)=[];
% error_mat=zeros(size(ox_data));
% error_mat(round(to(errors)))=1;
% to(errors)=[];
% errors=conv(error_mat,ones(1,1000),'same');
% noise_detector=nanmean(yo(fo>2 & fo<8,:),1)./nanmean(yo(fo>10 & fo<14,:),1);
% errors=conv(ox_data<=286 | ox_data==1023,ones(1,200)/200,'same')>0 | conv([diff(ox_data).^2;0],ones(1,200)/200,'same')>200;
% errors=conv(errors,ones(1000,1),'same')>0;

% ox_data(errors)=NaN;

% ox_data(isnan(ox_data))=interp1(find(~isnan(ox_data)),ox_data(~isnan(ox_data)),find(isnan(ox_data)),'linear');

% [buffer2]=BandFilt_Order(ox_data,fs,300,5,8);

[peaks, troughs,p_amp,t_amp]=findpeaks(buffer2(51:end));
p_amp_full=interp1(peaks,p_amp,1:length(ox_data),'linear','extrap');
p_amp_full(errors)=interp1(find(~errors),medfilt1(p_amp_full(~errors),1000),find(errors),'linear','extrap');
pa_mm=movmean(p_amp_full,1000);
lower_bound=pa_mm-2*nanstd(p_amp_full./pa_mm).*pa_mm;
binary_p=nan(1,length(ox_data));
binary_p(peaks)=p_amp;
binary_p(binary_p<lower_bound)=NaN;
peaks=find(~isnan(binary_p));
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
HR_full_interp=medfilt1(HR_full,5000,'omitnan');
HR_full2(isnan(HR_full))=interp1(find(~isnan(HR_full)),HR_full_interp(~isnan(HR_full)),find(isnan(HR_full)),'linear');
HR_full2(isnan(HR_full2))=interp1(find(~isnan(HR_full2)),HR_full2(~isnan(HR_full2)),find(isnan(HR_full2)),'nearest','extrap');
peaks_mat=binary_mat2;
% peaks_mat=peaks_mat./nanmean(peaks_mat);
HRV=sqrt(movmean(peaks_mat.^2,10000,'omitnan'));

HRV_interp=medfilt1(HRV,30000,'omitnan');
HRV(isnan(HRV))=interp1(find(~isnan(HRV)),HRV_interp(~isnan(HRV)),find(isnan(HRV)),'linear');
HRV(isnan(HRV))=interp1(find(~isnan(HRV)),HRV(~isnan(HRV)),find(isnan(HRV)),'nearest','extrap');
HRV=movmean(HRV,300);
data_struct.HR_fix=HR_full;
data_struct.HR_fix2=HR_full2;
data_struct.HRV=HRV;
data_struct.HRpeak_vector=peaks_mat;
% HR_full2=LowFilt_Order(HR_full,fs,30000,.5);
% plot(HR_full2)
% data_struct.HR_fix=HR_full2;