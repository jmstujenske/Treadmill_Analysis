% file='C:\Users\Admin\Dropbox (ListonLab)\ListonLab Team Folder\Joe\Animal2.tif';
% video_stack=bigread4(file);
function pupil_extract(path,correction,bg)
correction=correction-max(correction);
bg=bg-max(bg(:));
% path='H:\Treadmill_imaging_June2022\videos\grm2_mcherry_M4_FCRecall_d2.avi.mp4';
[directory,file,ext]=fileparts(path);
file=[file,ext];
% directory='H:\Treadmill_imaging_June2022\videos\';
% file='M2_day4.mp4';
filename=fullfile(directory,file);
% vr=vision.VideoFileReader(filename);
% matlab.video.read.UseHardwareAcceleration('on');
% reset(vr);
% video_frame=step(vr);
vr=VideoReader(filename);
% video_stack=zeros(size(video_frame,1),size(video_frame,2),4);
% video_stack(:,:,1)=video_frame(:,:,1);
% for a=2:4
%     video_frame=step(vr);
%     video_stack(:,:,a)=video_frame(:,:,1);
% end

video_stack=uint16(read(vr,[1 4]));video_stack=squeeze(video_stack(:,:,1,:))+uint16(reshape(-correction(1:4),[1 1 1 4]))+uint16(-bg);
input=pupil_calibration(nanmean(video_stack(:,:,1:4),3));
center=NaN;
radius=NaN;

% figure;
% reset(vr);
% for a=1:size(video_stack,3)
nF=vr.NumFrames;
data_save=zeros(nF,5);
figure;
for a=0:1000:nF
%     video_frame=step(vr);
video_frames=uint16(read(vr,[a+1 min(a+1000,nF)]))+uint16(reshape(-correction(a+1:min(a+1000,nF)),1,1,1,[]))+uint16(-bg);
% video_frames=squeeze(video_frames(:,:,1,:));
%     if ~isempty(video_frame)
%         count=count+1;
%     video_frame=uint8(video_frame*255);
parfor b=1:size(video_frames,4)
    [center,radius,area]=findpupil_inframe(video_frames(:,:,1,b),input);
    data_save(b+a,:)=[center radius area];
end
    imagesc(video_frames(input.eye_window(3):input.eye_window(4),input.eye_window(1):input.eye_window(2),1,end));
    hold on
    try
    rectangle('Position',[data_save(a+size(video_frames,4),1:2)-data_save(a+size(video_frames,4),3:4) data_save(a+size(video_frames,4),3:4)*2],'Curvature',1,'EdgeColor','r');
    end
    pause(0.02);
drawnow;
%     else
%         break
%     end
disp([num2str(min(a+1000,nF)),' Frames of ',num2str(nF),' Analyzed.'])
end
% data_save(count+1:end,:)=[];
% full_r=nanmean(data_save(:,[3 4]),2);
tointerp=false(1,size(data_save,1));
for a=1:4
    full_r=data_save(:,a);
limits=medfilt1(full_r,30);
tointerp=isnan(limits);
limits(tointerp)=interp1(find(~tointerp),limits(~tointerp),find(tointerp));
% blinks=logical(conv(full_r<limits | isnan(data_save(:,3)),ones(1,3),'same'));
blinks=((abs(full_r-limits))>(nanmean(abs(full_r-limits))+nanstd(abs(full_r-limits))*2)) | isnan(data_save(:,3));
blinks=logical(conv(blinks,ones(1,3),'same'));
% blinks=find(blinks);
tointerp=blinks |tointerp;
end
full_r=data_save(:,3);
% full_r(isnan(full_r))=0;
limits=movmedian(full_r,300,'omitnan')-4*movstd(full_r,300,'omitnan');
blinks=logical(conv(full_r<limits | isnan(data_save(:,3)),ones(1,3),'same'));
blinks=find(blinks);
data_save_new=data_save;
for a=1:4
data_save_new(tointerp,a)=interp1(find(~tointerp),data_save_new(~tointerp,a),find(tointerp));
end
data_save_new(:,5)=pi*prod(data_save_new(:,3:4),2);
radius=data_save_new(:,[3 4]);
center=data_save_new(:,[1 2]);
area=data_save_new(:,5);
save([directory,file(1:end-4),'_pupildata.mat'],'data_save','blinks','radius','center','area')
