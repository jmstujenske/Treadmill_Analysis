function data_out=extractdata_treadmill_timestamps_HR(timestamp_path)
if isstruct(timestamp_path)
    timestamp_path2=[];
    for a=1:length(timestamp_path)
        timestamp_path2{a}=dir2file(timestamp_path(a));
    end
else
    timestamp_path2={timestamp_path};
end
data_mat=[];
linecount=0;
for stamp_rep=1:length(timestamp_path2)
    fid=fopen(timestamp_path2{stamp_rep},'r');
    fseek(fid,0,'eof');
    flen=ftell(fid);
    fseek(fid,0,'bof');
    
    alldata=fread(fid,inf,'uint8');fclose(fid);
    first_data=find(alldata~=0,1,'first');
    first_to_read=first_data-7;
    if first_to_read<1
        first_to_read=1;
    end
    % first_to_read=1;
    numberoflines=floor((flen-first_to_read+1)/8);
    alldata=reshape(alldata(first_to_read:numberoflines*8+first_to_read-1),8,[]);
    data_mat=[data_mat;zeros(numberoflines,11)];
    if stamp_rep>1
        time1=datestr2num(timestamp_path(stamp_rep-1).date(end-7:end));
        time2=datestr2num(timestamp_path(stamp_rep).date(end-7:end));
        offset=round((time2-time1));
    else
        offset=0;
    end
    for a=1:numberoflines
        events=alldata(1,a);
        rot_value=alldata(2,a);
        mouseoxval=typecast(uint8(alldata([4 3],a)),'int16');
        timestamp=typecast(uint8(alldata([8:-1:5],a)),'uint32');
        todisplay=dec2bin(events);
        vecdisplay=zeros(1,8);
        vecdisplay(8-length(todisplay)+regexp(todisplay,'1','start'))=1;
        data_mat(a+linecount,1)=timestamp;
        data_mat(a+linecount,2)=rot_value;
        data_mat(a+linecount,3)=mouseoxval;
        data_mat(a+linecount,4:end)=vecdisplay;
    end
    
    %%sometimes the arduino messes up slightly. I am not sure why -- need to
    %%fix these "blips"
    
    ts=double(data_mat(max(linecount+1,1):linecount+numberoflines,1));
    in=1+find(diff(ts)<0 |  diff(ts)>2000,1,'first');
    if ~isempty(in)
    ts(in:end,1)=(ts(in-1)+1000*(1:(length(ts)-in+1))');
    end
    fix_timestamps_1=[false;diff(ts)>5000];
    fix_timestamps_2=[diff(ts)<=0;false];
    ends=find(fix_timestamps_2);
    starts=find(fix_timestamps_1);
    if ~isempty(ends)
    if ends(1)==1;
        starts=[1;starts];
    end
    end
    timestamps_fix=ts;
    if length(starts)==length(ends) && all(ends-starts)>=0
        fix_timestamps=[starts';ends'];
        for segments=1:size(fix_timestamps,2)
            timestamps_fix(fix_timestamps(1,segments)-1:fix_timestamps(2,segments)+1)=...
                linspace(ts(fix_timestamps(1,segments)-1),ts(fix_timestamps(2,segments)+1),fix_timestamps(2,segments)-fix_timestamps(1,segments)+3);
        end
    end
    data_mat(max(linecount+1,1):linecount+numberoflines,1)=timestamps_fix;
    linecount=linecount+numberoflines;
    %%timestamps fixed
    if stamp_rep>1
        
        data_mat(end-numberoflines+1:end,1)=data_mat(end-numberoflines+1:end,1)-data_mat(end-numberoflines+1,1)+data_mat(1,1)+offset*1e6;
    end
end
% try
% ox_data=data_mat(:,3);
% ox_data=sgolayfilt(ox_data,3,21);
% fs=1000;
% [yo, fo, to]=mtcsg(ox_data./sqrt(movmean(ox_data.^2,20000)),fs*5,fs,fs,fs-500,1.5);
% [~,in]=max(yo(fo>2 &fo<=15.5,:).*fo(fo>2 &fo<=15.5));
% to=to*fs+fs/2;
% 
% error_index=(movmean(nanmean(yo(fo>2 & fo<5,:).*fo(fo>2 & fo<5),1),1));
% error_index=error_index-movmin(error_index,50);
% 
% error_index=sqrt(error_index);
% errors=error_index>1.5*nanstd(error_index);
% errors=errors' | fo(in+find(fo<=2,1,'last'))<5;
% 
% error_mat=zeros(size(ox_data));
% error_mat(round(to(errors)))=1;
% 
% errors=conv(error_mat,ones(1,1000),'same');
% 
% [buffer2]=BandFilt_Order(ox_data,fs,100,5,15);
% [peaks, troughs,p_amp,t_amp]=findpeaks(buffer2(51:end));
% peaks=peaks(p_amp>=5);
% peaks=peaks+50;
% peaks(find(diff(peaks)<10)+1)=[];
% % peaks(ismember(peaks,find(errors)))=NaN;
% HR=1e6./[double(data_mat(peaks(5:end),1))-double(data_mat(peaks(1:end-4),1))]*60*4;
% HR(HR<300 | HR>900)=NaN;
% tofix=ismember(peaks(3:end-1),find(errors));
% HR(tofix)=NaN;
% tofix=isnan(HR);
% peaktimes=peaks(3:end-2);
% HR_full=interp1(peaktimes(~tofix),HR(~tofix),1:length(ox_data),'linear');
% peaktimes(tofix)=NaN;
% tofix=abs(HR_full-medfilt1(HR_full,5000))>20;
% peaktimes(ismember(peaktimes,find(tofix)))=NaN;
% HR_full=interp1(find(~tofix),HR_full(~tofix),1:length(ox_data),'linear');
% HR_full(isnan(HR_full))=interp1(find(~isnan(HR_full)),HR_full(~isnan(HR_full)),find(isnan(HR_full)),'linear');
% HR_full2=LowFilt_Order(HR_full,fs,30000,.5);
% data_out.HR=HR_full2;
% data_out.pulse_ox_peaks=peaktimes;
% data_out.HR=uint16(data_out.HR);
% catch
%     data_out.HR=[];
%     data_out.pulse_ox_peaks=[];
% end
data_out.timestamps=data_mat(:,1)-data_mat(1,1);

%%fix when clock resets
sample_diff=diff((data_out.timestamps));
discont=find(sample_diff<-1e9)+1;
if ~isempty(discont)
samp_space=(median(sample_diff));
for dis_rep=1:length(discont)
    data_out.timestamps(discont(dis_rep):end)=data_out.timestamps(discont(dis_rep):end)-data_out.timestamps(discont(dis_rep))+data_out.timestamps(discont(dis_rep)-1)+samp_space;
end
end
data_out.position=data_mat(:,2);
data_out.pulse_ox_raw=data_mat(:,3);
data_out.position=uint8(data_out.position);
data_out.timestamps=uint32(data_out.timestamps);
data_out.pulse_ox_raw=uint16(data_out.pulse_ox_raw);

thermal_cam_in=find(all(logical(data_mat),2));
data_out.thermal_cam_ontime=data_out.timestamps(thermal_cam_in);
data_mat(thermal_cam_in,4:end)=0;
data_out.twophoton_frame_times=data_out.timestamps(data_mat(:,5)>0);
licks=zeros(size(data_mat(:,6)));
licks(logical(data_mat(:,6)))=1;
licks(logical(data_mat(:,7)))=2;
licks(licks==0)=[];
data_out.lick_type=licks;
data_out.sounds_in=find(data_mat(:,8));
data_out.random_reward_times=data_out.timestamps(data_mat(:,9)>0);
data_out.shocks_in=find(data_mat(:,10));
data_out.video_frame_times=data_out.timestamps(data_mat(:,4)>0);
data_out.licks_times=data_out.timestamps(data_mat(:,6) | data_mat(:,7));
dts=diff(data_out.video_frame_times);
data_out.video_sample_rate=1e6./nanmean(dts(dts<40000 | dts>20000));
data_out.event_times=data_out.timestamps(data_mat(:,8)>0);
data_out.shock_times=data_out.timestamps(data_mat(:,10)>0);
data_out.movement_times=data_out.timestamps(data_mat(:,11)>0);