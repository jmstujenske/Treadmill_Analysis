function data_out=extract_pupildata(data_out);

fid=fopen(data_out.bin_path,'r');
pupil=fread(fid,inf,'uint16');
fclose(fid);
if length(data_out.video_frame_times)<20
    vrobj=VideoReader([data_out.bin_path(1:end-10),'.avi.mp4']);
    nf=vrobj.NumFrames;
    lastframe=double(data_out.timestamps(end))-110000;
    ivi=33295;
    frametimes=(0:nf-1)*ivi;
    frametimes=frametimes-frametimes(end)+lastframe;
    data_out.video_frame_times=uint32(frametimes)';
%     keyboard
end
space=mode(diff(find(diff(pupil)./pupil(1:end-1)>10)));
remainder=mod(size(pupil,1),space);
lastbit=pupil(end-remainder+1:end);
pupil=pupil(1:end-remainder);
pupil=reshape(pupil,12,5,[]);
remainder2=mod(size(lastbit,1),5);
lastbit=lastbit(1:end-remainder2);
pupil_lastbit=reshape(lastbit,[],5);
pupil2=reshape(permute(pupil,[1 3 2]),size(pupil,1)*size(pupil,3),5);
pupil2=cat(1,pupil2,pupil_lastbit);
pupil2(pupil2==0)=NaN;
if length(data_out.video_frame_times)~=length(pupil2)
    frame_times_new=zeros(1,length(pupil2));
    find_breaks=find(diff(data_out.video_frame_times)>40000);
    dts=diff(data_out.video_frame_times);
    acq_rate=nanmean(dts(dts<40000 & dts>30000));
    for break_rep=1:length(find_breaks)
        num2add=double(ceil((data_out.video_frame_times(find_breaks(break_rep)+1)-data_out.video_frame_times(find_breaks(break_rep)))/acq_rate));
        num2add=num2add+1;
        data_out.video_frame_times=[data_out.video_frame_times(1:find_breaks(break_rep)-1);...
            uint32(linspace(double(data_out.video_frame_times(find_breaks(break_rep))),...
            double(data_out.video_frame_times(find_breaks(break_rep)+1)),num2add))';...
            data_out.video_frame_times(find_breaks(break_rep)+2:end)];
        find_breaks(break_rep+1:end)=find_breaks(break_rep+1:end)+num2add-2;
    end
    
end
data_out.pupil.times=data_out.video_frame_times(1:min(length(pupil2),length(data_out.video_frame_times)));
% pupil2=pupil2(1:min(length(pupil2),length(data_out.video_frame_times)),:);

pupil_center=nanmedian(pupil2(:,[1 2]));
pupil_radius=nanmedian(pupil2(:,[3 4]));
% pupil_radius_upper=pupil_radius+2.5*nanstd(pupil_radius);

% upper_lim=pupil_center+pupil_radius_upper;
% lower_lim=pupil_center-pupil_radius_upper;
% pupil2(pupil2(:,1)>upper_lim(1),:)=NaN;
% pupil2(pupil2(:,2)>upper_lim(2),:)=NaN;
% nanmean(data_out.pupil.area)-2.5*nanstd(data_out.pupil.area)
data_out.pupil.area=pupil2(:,5);
lower=nanmean(data_out.pupil.area)-3*nanstd(data_out.pupil.area);
upper=nanmean(data_out.pupil.area)+5*nanstd(data_out.pupil.area);
data_out.pupil.area(data_out.pupil.area<lower | data_out.pupil.area>upper)=NaN;
while 1
    ins=find(~isnan(data_out.pupil.area));
    if ~isempty(data_out.pupil.times)
    offsets=diff(data_out.pupil.area(ins))./diff(double(data_out.pupil.times(ins)));
    else
        offsets=diff(data_out.pupil.area(ins))./(diff(ins)*33000);
    end
    toremove=ins(find(abs(offsets)>.004)+1);
    if ~isempty(toremove)
    data_out.pupil.area(toremove)=NaN;
    else
        break;
    end
end
pupil2(isnan(data_out.pupil.area),:)=NaN;
data_out.pupil.coords=pupil2(:,[1 2]);
data_out.pupil.radius=pupil2(:,[3 4]);
tofix=isnan(data_out.pupil.area) | any(data_out.pupil.radius>200,2);
data_out.pupil.area(tofix)=interp1(find(~tofix),data_out.pupil.area(~tofix),find(tofix));
data_out.pupil.coords(tofix,1)=interp1(find(~tofix),data_out.pupil.coords(~tofix,1),find(tofix));
data_out.pupil.coords(tofix,2)=interp1(find(~tofix),data_out.pupil.coords(~tofix,2),find(tofix));
data_out.pupil.radius(tofix,1)=interp1(find(~tofix),data_out.pupil.radius(~tofix,1),find(tofix));
data_out.pupil.radius(tofix,2)=interp1(find(~tofix),data_out.pupil.radius(~tofix,2),find(tofix));
data_out.pupil.area=uint16(data_out.pupil.area);
data_out.pupil.coords=uint16(data_out.pupil.coords*100);
data_out.pupil.radius=uint16(data_out.pupil.radius*100);