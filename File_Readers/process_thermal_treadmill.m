function thermal_data=process_thermal_treadmill(filename,cam_start_time,frames)

if nargin<2 || isempty(cam_start_time)
    cam_start_time=0;
end
if nargin<3 || isempty(frames)
    frames=[1 inf];
end
if length(frames)==1
    if ~isinf(frames(1))
    frames=repmat(frames,1,2);
    else
        frames=[1 inf];
    end
end
fr=28;
subsample_factor=5;
init_batch=500;
if diff(frames)+1<fr*2
    error('Not enough frames.');
end

if diff(frames)+1<init_batch*4
    staggered=false;
else
    staggered=true;
end
% filename='H:\Treadmill_imaging_June2022\videos\grm2_mcherry_M1_licktrain3_thermal';


if isstr(filename)
[path,name,ext]=fileparts(filename);
path=[path,'\'];
filename=name;
bin_file=[path,filename,'.bin'];
frame_file=[path,filename,'.txt'];
else
    bin_in=cellfun(@(x)contains(x,'.bin'),{filename(:).name});
    frame_in=cellfun(@(x)contains(x,'.txt'),{filename(:).name});
    bin_file=fullfile(filename(bin_in).folder,filename(bin_in).name);
    frame_file=fullfile(filename(frame_in).folder,filename(frame_in).name);
end
fid=fopen(bin_file,'r');
[path,name,ext]=fileparts(bin_file);
fseek(fid,0,'eof');
flen=ftell(fid);
fclose(fid);
f_w=320;
f_h=240;
nf_tot=flen/f_w/f_h/2;
frames(isinf(frames))=nf_tot;
nf=diff(frames)+1;
    frame_stamps=importdata(frame_file);
frame_stamps(:,1)=frame_stamps(:,1)+1;
% frame_stamps(:,1)=1:length(frame_stamps(:,1));
frame_stamps(:,2)=(frame_stamps(:,2)-frame_stamps(1,2))*1e3+double(cam_start_time);
frame_stamps=frame_stamps(frames(1):frames(2),:);
if ~exist(fullfile(path,[name,'_temp.bin']))

m_fid=memmapfile(bin_file,'Format',{'uint16' [f_w f_h nf_tot] 'data'},'Writable',true);
firstn=1000;
tofixfirst=double(m_fid.Data.data(:,:,min(1:firstn,nf_tot)))-imgaussfilt(double(m_fid.Data.data(:,:,min(1:firstn,nf_tot))),10)<-10000;
tocorrect=m_fid.Data.data(:,:,min(1:firstn,nf_tot))<10000 & tofixfirst;
if any(tocorrect,1:3)
data_sub=m_fid.Data.data(:,:,min(1:firstn,nf_tot));
data_sub(tocorrect)=2^16-1+data_sub(tocorrect);
m_fid.Data.data(:,:,min(1:5000,nf_tot))=data_sub;
end
checkfix=find(m_fid.Data.data(:,:,1001:end)<10000);
if ~isempty(checkfix)
    m_fid.Data.data(checkfix)=2^16-1+m_fid.Data.data(checkfix);
end

data_downsample=imresize(m_fid.Data.data(:,:,frames(1):frames(2)),1/subsample_factor);
% options_nonrigid = NoRMCorreSetParms('d1',size(data_downsample,1),'d2',size(data_downsample,2),'grid_size',[12 12],'overlap_pre',[6 6],'mot_uf',4,'bin_width',3000,'max_shift',.2,'max_dev',[20 20],'us_fac',50,'init_batch',500);

options_rigid = NoRMCorreSetParms('d1',size(data_downsample,1),'d2',size(data_downsample,2),'bin_width',5000,'max_shift',10,'us_fac',50,'init_batch',init_batch,'print_msg',false);
% nf=size(data_downsample,3);

trend=quantile(reshape(data_downsample,[],size(data_downsample,3)),.5,1);
trend=trend./nanmean(trend);
trend(trend<.01)=.01;
data_downsample=data_downsample./reshape(trend,1,1,[]);
% cleaned_stack=((data_downsample)-imclose(data_downsample,strel('disk',15)));
% data_ds=gpuArray(data_downsample);
block_size=10000;
nblock=nf/block_size;
cleaned_stack=zeros(size(data_downsample),class(data_downsample));
for block_rep=1:nblock
data_ds=gpuArray(uint8(data_downsample(:,:,min((1:block_size)+(block_rep-1)*block_size,nf))/(2^8)));
cleaned_stack(:,:,min((1:block_size)+(block_rep-1)*block_size,nf))=double(gather(imopen(uint8(double(data_ds)-double(imclose(data_ds,strel('disk',20)))+100),strel('disk',1))));
end
% cleaned_stack=imopen((data_downsample)-imclose(data_downsample,strel('disk',20)),strel('disk',1));
% cleaned_stack=((data_downsample)-imclose(data_downsample,strel('disk',20)));
init_batch=min(init_batch,nf);
if staggered
tic; [~,shifts1,template1,options_rigid] = normcorre(single(cat(3,cleaned_stack(:,:,init_batch+1:end),cleaned_stack(:,:,1:init_batch*2))),options_rigid); toc
% motcorr_stack=cat(3,motcorr_stack(:,:,end-init_batch*2+1:end),motcorr_stack(:,:,init_batch+1:end-init_batch*2));
shifts1=shifts1([end-init_batch*2+1:end init_batch+1:end-init_batch*2]);
else
tic; [~,shifts1,template1,options_rigid] = normcorre(single(cleaned_stack),options_rigid); toc
end
options_rigid = NoRMCorreSetParms('d1',size(m_fid.Data.data,1),'d2',size(m_fid.Data.data,2),'bin_width',5000,'max_shift',10,'us_fac',50,'init_batch',500,'print_msg',false);

for a=1:length(shifts1)
    f_names={'shifts','shifts_up','diff'};
    for b=1:3
    shifts1(a).(f_names{b})=shifts1(a).(f_names{b})*subsample_factor;
%     shifts2(a).(f_names{b})=shifts2(a).(f_names{b})*subsample_factor;
    end
end
chunk_size=20000;
fid_w=fopen(fullfile(path,[name,'_temp.bin']),'w');
n_chunks=ceil((diff(frames)+1)/chunk_size);
for chunk_rep=1:n_chunks
block_size=1000;
if chunk_rep==n_chunks
    remainder=mod(diff(frames)+1,chunk_size);
    remainder(remainder==0)=chunk_size;
else
    remainder=chunk_size;
end
M_final=zeros([size(m_fid.Data.data,1:2),remainder],'single');
n_blocks=ceil(remainder/block_size);
% options_nonrigid = NoRMCorreSetParms('d1',size(data,1),'d2',size(data,2),'grid_size',[12 12]*5,'overlap_pre',[6 6]*5,'mot_uf',4,'bin_width',3000,'max_shift',.2,'max_dev',[20 20],'us_fac',50,'init_batch',500);

for a=1:n_blocks
    M_final(:,:,min((a-1)*block_size+(1:1000),remainder))=apply_shifts(m_fid.Data.data(:,:,min((a-1)*block_size+(1:1000)+frames(1)-1+(chunk_rep-1)*chunk_size,nf_tot)),shifts1(min((a-1)*block_size+(1:1000)+(chunk_rep-1)*chunk_size,nf)),options_rigid);
end
% M_final=motcorr_stack;
M_final=permute(M_final,[2 1 3]);
M_final=M_final(subsample_factor*3+1:end-subsample_factor*3,subsample_factor*3+1:end-subsample_factor*3,:);
fwrite(fid_w,M_final,'single');
end
fclose(fid_w);
end
M_final_fid=memmapfile(fullfile(path,[name,'_temp.bin']),'Format',{'single' [f_h-subsample_factor*6 f_w-subsample_factor*6 nf] 'data'},'Writable',true);
clear M_final
trend=(squeeze(nanmean(M_final_fid.Data.data,1:2)));
trend=medfilt1(trend,fr,'truncate');
M_final_fid.Data.data=(M_final_fid.Data.data./reshape(trend,1,1,[])).*nanmean(trend);
% M_final=double(M_final);

% mask=nanstd(gpuArray(dd),[],3);
startin=max(5000-frames(1)+1,1);
subin=5;
if startin*2>size(M_final_fid.Data.data,3)
    startin=1;
    subin=1;
end
dd=M_final_fid.Data.data(:,:,(startin:subin:end-subin)+1)-M_final_fid.Data.data(:,:,startin:subin:end-subin);
mask=nanstd(dd,[],3);
mask(:,[1:subsample_factor*5 end-subsample_factor*5+1:end])=NaN;
mask([1:subsample_factor*3 end-subsample_factor*3+1:end],:)=NaN;

mask=(mask-max(nanmin(mask(:)),0))./(nanmax(mask(:))-max(nanmin(mask(:)),0));
% mask2=mask;
% mask2=(mask2-min(mask2(:)))./(max(mask2(:))-min(mask2(:)));
mask=mask>.5;
dd=diff(M_final_fid.Data.data(:,:,startin:min(startin+chunk_size,nf)),[],3);
trace=squeeze(nanmean(dd.*mask,1:2))/sum(mask(:))./nanmean(mask(:));
trace=trace-movmedian(trace,fr);
% trace=gather_try(trace);
tofix=trace>10;
trace(tofix)=interp1(find(~tofix),trace(~tofix),find(tofix),'linear','extrap');
A=corr(reshape(dd,[],size(dd,3))',BandFilt_Order(trace,fr,fr,1,6));
mask=reshape(A,size(M_final_fid.Data.data,1:2));
mask((mask)<max(A(:))*.5)=NaN;
mask_bg=reshape(A,size(M_final_fid.Data.data,1:2));
mask_bg((mask_bg)>max(A(:))*.3)=NaN;

mask(isnan(mask))=0;
% mask=gather_try(mask);
% mask_bg=gather_try(mask_bg);
valid=mask>0;
valid(:,[1:subsample_factor*8 end-subsample_factor*8+1:end])=0;
CC=bwconncomp(valid);
mask_orig=mask;
mask_orig=mask_orig>0;
if length(CC.PixelIdxList)>1
patch_sizes=cellfun(@numel,CC.PixelIdxList);
        [~,size_in]=sort(patch_sizes,'descend');
largest2=size_in(1:2);
valid_mask=zeros(size(mask_orig));
for a=1:2
    valid_mask(CC.PixelIdxList{largest2(a)})=a;
end
mask=mask_orig.*(valid_mask>0);
center=regionprops(CC,'Centroid');
center=center(largest2);
if center(1).Centroid(1)<center(2).Centroid(1)
mask_l=mask_orig.*(valid_mask==1);
mask_r=mask_orig.*(valid_mask==2);
else
    mask_l=mask_orig.*(valid_mask==2);
    mask_r=mask_orig.*(valid_mask==1);
end
else
    mask_l=mask_orig;
    mask_r=zeros(size(mask_orig));
end
% mask=mask>0;
mask_full=mask./sum(mask(:));
mask_l=mask_l/sum(mask_l(:));
mask_r=mask_r/sum(mask_r(:));
mask_bg=~isnan(mask_bg);
mask_bg=mask_bg./sum(mask_bg(:));
thermal_data=struct;
%for some reason, python sometimes spits out extraneous timestamps that do
%not correspond to anything. This is indicated by repeated frame numbers
%The extraneous timestamps seem to come out at the expected frequency of
%acquisition, which is odd. I don't know the cause of this, but I have
%written a fix. May need to change the shutter_interval, if it is changed
%in the future. Right now, it is set to calibrate every 400 frames after
%the beginning.
shutter_interval=400;
ins_shutter=find(diff(frame_stamps(:,2))>100e3);
which_haveshutter=mode(mod(ins_shutter,shutter_interval));
[~,in]=unique(frame_stamps(:,1));
frame_stamps_fixed=frame_stamps(in,:);
intervals=diff(frame_stamps_fixed(:,2));
typical_interval=nanmean(intervals(intervals<40e3));
typical_shutter=nanmean(intervals(intervals>200e3 & intervals<1000e3)); %will fix timestamps with what is expected for 
breaks=find(diff(frame_stamps_fixed(:,2))>1000e3);
for b=1:length(breaks)
    time_diff=diff(frame_stamps_fixed([breaks(b) breaks(b)+1],2));
    hasshutter=mod(frame_stamps_fixed(breaks(b),1),shutter_interval)==which_haveshutter;
    if ~hasshutter
        frame_stamps_fixed(breaks(b)+1:end,2)=frame_stamps_fixed(breaks(b)+1:end,2)-time_diff+typical_interval;
    else
                frame_stamps_fixed(breaks(b)+1:end,2)=frame_stamps_fixed(breaks(b)+1:end,2)-time_diff+typical_shutter;
    end
end
frame_stamps=frame_stamps_fixed;
%%

for mask_rep=1:3
    switch mask_rep
        case 1
            mask=mask_full;
        case 2
            mask=mask_l;
        case 3
            mask=mask_r;
    end

M_final=sparse([]);
n_chunks=ceil(nf/chunk_size);
for chunk_rep=1:n_chunks
if chunk_rep==n_chunks
    remainder=mod(diff(frames)+1,chunk_size);
    remainder(remainder==0)=chunk_size;
else
    remainder=chunk_size;
end
M_final=cat(2,M_final,sparse(reshape(double(M_final_fid.Data.data(:,:,(1:remainder)+(chunk_rep-1)*chunk_size)).*mask,[],remainder)));
end
trace_alt=full(squeeze(nanmean(M_final,1)))';
% dd=diff(M_final.*mask,[],3);
subsample_time=10;
dd2=diff(M_final,[],2);

threshold=quantile(dd2(:,1:subsample_time:end),.15,2);
threshold2=quantile(dd2(:,1:subsample_time:end),.85,2);

trace2=nanmean(((dd2)<threshold)-((dd2)>threshold2),1);
trace2=full(squeeze(trace2)');

trace_alt=trace_alt-movmean(trace_alt,fr);


resample_factor=5;
trace_smooth=interp1(1:length(trace2),trace2,1:1/resample_factor:length(frame_stamps)+(resample_factor-1)/resample_factor,'pchip')';
alt_trace_smooth=interp1(1:length(trace_alt),trace_alt,1:1/resample_factor:length(frame_stamps)+(resample_factor-1)/resample_factor,'pchip')';

frame_upsample=interp1(frame_stamps(1:end,1),frame_stamps(1:end,2),1:1/resample_factor:length(frame_stamps)+(resample_factor-1)/resample_factor);
frame_upsample=[(1:1/resample_factor:length(frame_stamps)+(resample_factor-1)/resample_factor)' frame_upsample'];
movements=squeeze(nanmean(abs(dd.*mask_bg),1:2));
movements_smooth=interp1(1:length(movements),movements,1:1/resample_factor:length(frame_stamps)+(resample_factor-1)/resample_factor,'pchip')';

movements_smooth=LowFilt_Order(movements_smooth,fr*5,fr*5*3,.1);
% movements_smooth=[movements_smooth;zeros(resample_factor,1)];

toremove=isnan(frame_upsample(:,2));
trace_smooth(toremove)=[];
alt_trace_smooth(toremove)=[];
frame_upsample(toremove,:)=[];
movements_smooth(toremove,:)=[];
discont_start=find(diff(frame_stamps(:,2))>100*1e3);
% discont_stop=find(diff(frame_stamps(:,2))>100)+1;
for a=1:length(discont_start)
    [in1]=find(frame_upsample(:,1)==discont_start(a));
    [in2]=find(frame_upsample(:,1)==discont_start(a)+1);
    in1=in1+1;in2=in2-1;
    trace_smooth(in1:in2)=NaN;
    movements_smooth(in1:in2)=NaN;
    alt_trace_smooth(in1:in2)=NaN;
end
fps_new=100;
space=1e6/fps_new;
ts_even=frame_upsample(1,2):space:frame_upsample(end,2);
trace_smooth_even=interp1(frame_upsample(:,2),trace_smooth,ts_even);
alt_trace_smooth_even=interp1(frame_upsample(:,2),alt_trace_smooth,ts_even);
movements_even=interp1(frame_upsample(:,2),movements_smooth,ts_even);
% thermal_data(mask_rep)=struct;
thermal_data(mask_rep).mask=mask;
thermal_data(mask_rep).mask_bg=mask_bg;
thermal_data(mask_rep).inhale=trace_smooth_even';
thermal_data(mask_rep).temp=alt_trace_smooth_even';
thermal_data(mask_rep).ts=ts_even';
thermal_data(mask_rep).raw.ts=frame_stamps(:,2);
thermal_data(mask_rep).raw.inhale=[trace2;0];
thermal_data(mask_rep).raw.temp=[trace_alt];
thermal_data(mask_rep).raw.residual=movements;
thermal_data(mask_rep).residual=movements_even;
% thermal_data(mask_rep).temp_interp=thermal_data(mask_rep).temp;

% thermal_data(mask_rep).temp_interp(isnan(thermal_data(mask_rep).temp))=interp1(find(~isnan(thermal_data(mask_rep).temp)),...
%     thermal_data(mask_rep).temp(~isnan(thermal_data(mask_rep).temp)),find(isnan(thermal_data(mask_rep).temp)),'spline');
thermal_data(mask_rep).temp_interp=interp_periodic(thermal_data(mask_rep).temp,fps_new,[1 inf],4);
% thermal_data(mask_rep).inhale_interp=thermal_data(mask_rep).inhalations;
% thermal_data(mask_rep).inhale_interp(isnan(thermal_data(mask_rep).inhalations))=interp1(find(~isnan(thermal_data(mask_rep).inhalations)),thermal_data(mask_rep).inhalations(~isnan(thermal_data(mask_rep).inhalations)),...
%     find(isnan(thermal_data(mask_rep).inhalations)),'spline');
thermal_data(mask_rep).inhale_interp=interp_periodic(thermal_data(mask_rep).inhale,fps_new,[1 inf],4);

ins=find(histc(thermal_data(mask_rep).raw.ts,thermal_data(mask_rep).ts));
ins=[ins;length(thermal_data(mask_rep).ts)];

breaks=find(diff(ins)>20);
ts=thermal_data(mask_rep).raw.ts;
ts_new=ts;
% raw.temp_interp=thermal_data(mask_rep).raw.temp;
% raw.inhale_intepr=thermal_data(mask_rep).raw.inhale;
averagediff=3.78;
addedsamples=0;
for rep=1:length(breaks)
    curr_break=breaks(rep)+addedsamples;
    break_length=diff(ins([breaks(rep) breaks(rep)+1]));
    n_toadd=ceil(break_length/averagediff)-1;
    toadd=linspace(ts(breaks(rep)),ts(breaks(rep)+1),n_toadd+2);
    toadd=toadd(2:end-1)';
    addedsamples=addedsamples+n_toadd;
    ts_new=[ts_new(1:curr_break);toadd;ts_new(curr_break+1:end)];
    
end

ins_new=find(histc(ts_new,thermal_data(mask_rep).ts));
ins_new=[ins_new;length(thermal_data(mask_rep).ts)];
thermal_data(mask_rep).raw.temp_interp=thermal_data(mask_rep).temp_interp(ins_new);
thermal_data(mask_rep).raw.inhale_interp=thermal_data(mask_rep).inhale_interp(ins_new);
thermal_data(mask_rep).raw.ts_interp=ts_new;
starts=find(diff(thermal_data(mask_rep).raw.inhale_interp>0)==1)+1;
ends=find(diff(thermal_data(mask_rep).raw.inhale_interp>0)==-1);

if starts(1)>ends(1)
    ends(1)=[];
end
if starts(end)>ends(end)
    starts(end)=[];
end
thermal_data(mask_rep).raw.inhale_stats.starts_idx=starts;
thermal_data(mask_rep).inhale_stats.starts_idx=ins_new(starts);
thermal_data(mask_rep).raw.inhale_stats.ends_idx=ends;
thermal_data(mask_rep).inhale_stats.ends_idx=ins_new(ends);
n=length(thermal_data(mask_rep).inhale_stats.starts_idx);
% thermal_data(mask_rep).inhale_size=zeros(1,n);
thermal_data(mask_rep).inhale_stats.AUC=zeros(1,n);
thermal_data(mask_rep).inhale_stats.tempchange=zeros(1,n);
thermal_data(mask_rep).inhale_stats.length=zeros(1,n);
for inhale_rep=1:n
    val=max(thermal_data(mask_rep).raw.inhale_interp(thermal_data(mask_rep).raw.inhale_stats.starts_idx(inhale_rep):thermal_data(mask_rep).raw.inhale_stats.ends_idx(inhale_rep)));
    thermal_data(mask_rep).inhale_peak(inhale_rep)=val;
        val=sum(thermal_data(mask_rep).raw.inhale_interp(thermal_data(mask_rep).raw.inhale_stats.starts_idx(inhale_rep):thermal_data(mask_rep).raw.inhale_stats.ends_idx(inhale_rep)));
    thermal_data(mask_rep).inhale_stats.AUC(inhale_rep)=val;
    if inhale_rep~=1
    val1=max(thermal_data(mask_rep).raw.temp_interp([thermal_data(mask_rep).raw.inhale_stats.ends_idx(inhale_rep-1):thermal_data(mask_rep).raw.inhale_stats.ends_idx(inhale_rep)]));
    else
        val1=max(thermal_data(mask_rep).raw.temp_interp([1:thermal_data(mask_rep).raw.inhale_stats.ends_idx(inhale_rep)]));
    end
    if inhale_rep~=n
        val2=min(thermal_data(mask_rep).raw.temp_interp([thermal_data(mask_rep).raw.inhale_stats.starts_idx(inhale_rep):thermal_data(mask_rep).raw.inhale_stats.starts_idx(inhale_rep+1)]));
    else
       val2=min(thermal_data(mask_rep).raw.temp_interp(thermal_data(mask_rep).raw.inhale_stats.starts_idx(inhale_rep):end));
    end
    val=diff([val1 val2]);
    thermal_data(mask_rep).inhale_stats.tempchange(inhale_rep)=val;
        val=diff(thermal_data(mask_rep).ts([thermal_data(mask_rep).inhale_stats.starts_idx(inhale_rep),thermal_data(mask_rep).inhale_stats.ends_idx(inhale_rep)]));
    thermal_data(mask_rep).inhale_stats.length(inhale_rep)=val;
end
end

% save([path,filename,'.mat'],'thermal_data(mask_rep)');

% thermal_data(mask_rep).inhalations(thermal_data(mask_rep).m_corr_residual>500)=0;
% starts=find(diff(max(thermal_data(mask_rep).inhale_interp,0)>0)==1)+1;
% ends=find(diff(max(thermal_data(mask_rep).inhale_interp,0)>0)==-1);
% IBI=[NaN;starts(3:end)-starts(1:end-2);NaN];
% 
% peak_amp=zeros(1,length(starts));
% for a=1:length(starts)-1
%     peak_amp(a)=range(thermal_data(mask_rep).temperature_nostril_interp(starts(a):(ends(a)+starts(a+1))/2));
% end
% 
% peak_vel=zeros(1,length(starts));
% for a=1:length(starts)-1
%     peak_vel(a)=max(thermal_data(mask_rep).inhale_interp(starts(a):ends(a)));
% end
% 
% % bg=reshape(A,size(M_final,1:2)-4);
% % bg(bg>.2 |bg<0)=NaN;
% bg=mask2<1.5;
% % bwconncomp(mask2)
% trace=squeeze(nanmean(nanmean(diff(M_final(3:end-2,10:end-9,:).*mask,[],3),1),2))/sum(mask(:));
% bg_trace=squeeze(nanmean(nanmean(diff(M_final(3:end-2,10:end-9,:).*(bg),[],3),1),2))/sum(bg(:));
% trace_smooth=sgolayfilt(interp1(1:length(trace),double(trace),1:1/5:length(trace)),3,11);
% % bg_smooth=sgolayfilt(interp1(1:length(bg_trace),double(bg_trace),1:1/5:length(bg_trace)),3,21);
% % [U,S,V]=svdecon(reshape(diff(M_final,[],3),prod(size(M_final,1:2)),[]));
% % trace_smooth=sgolayfilt(interp1(1:length(trace),double(V(:,1)),1:1/5:length(trace)),3,21);
% [p,t,p_s,t_s]=findpeaks(trace_smooth);
% p=round((p-1)/5+1);
%     count=33100;
%     rr=[min(M_final(:))*.95 max(M_final(:))*.95];
%     figure;for a=thermal_data(mask_rep).ts_100Hz_thermal(count+1:end-1000)'
% 
%     count=count+1;
%         [~,in]=min(abs(thermal_data(mask_rep).raw.ts-a));
%         subplot(1,2,1);
%         imagesc(M_final(:,:,in));
% %         imagesc(diff(motcorr_stack(:,:,[in in+1]),[],3))
%         caxis(rr);
% %         caxis([-3000 3000]);
% %         caxis([0 400]);
%         subplot(1,2,2);
%         plot(1:201,thermal_data(mask_rep).temperature_nostril_interp(count-100:count+100));
%         hold on;
%         plot([101 101],[-.6 .6],'k-');
%         hold off;
%         drawnow;
%         pause(.001);
%     end

% [f1,f2] = freqspace(42,'meshgrid');
%         Hd = ones(42);
%         r = sqrt(f1.^2 + f2.^2);
%         Hd((r>.1)) = 0;
%         h = fwind1(Hd,hamming(42));
%         freqz2(h);
 