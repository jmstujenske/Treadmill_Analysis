function input=pupil_calibration(vid_file)
TimeStamp_directory='C:\FearConditioningData\TimeStamps\';
Video_directory='C:\FearConditioningData\Videos\';

% [TS_file,TS_path]=uigetfile([TimeStamp_directory,'*.txt'], 'Select the TimeStamp file');
if nargin<1 || isempty(vid_file)
[vid_file,vid_path]=uigetfile([Video_directory,'*.avi'], 'Select the Video file');
vid_file=[vid_path,vid_file];
end
if ischar(vid_file)
vidObj=VideoReader(vid_file);
template=read(vidObj,[1 30]);
else
    template=vid_file;
end

template=nanmean(template(:,:,1,:),4);
intensityrange=double([0 max(template(:))]);
f=uifigure('Position',[100 100 700 700],'Resize','off');
% p=uipanel(f,'Position',[30 300 500 300]);
ax=uiaxes(f,'Position',[30 300 480 360]);
default_th=ceil(30/255*intensityrange(2));
slider = uispinner(f,'Position',[80 295 100 22],'Limits',intensityrange,'Value',default_th,'ValueChangedFcn',@(sld,event) updateColorMap(sld,[],ax,intensityrange),'ValueChangingFcn',@(sld,event) updateColorMap(sld,[],ax,intensityrange));
slider2 = uispinner(f,'Position',[200 295 100 22],'Limits',[-180 180],'Value',0,'ValueChangedFcn',@(sld,event) sld_rotate(sld,slider,ax,intensityrange),'ValueChangingFcn',@(sld,event) sld_rotate(sld,slider,ax));
set(slider,'ValueChangedFcn',@(sld,event) updateColorMap(sld,slider2,ax,intensityrange),'ValueChangingFcn',@(sld,event) updateColorMap(sld,slider2,ax,intensityrange));

button = uibutton(f,'Position',[80 255 100 30],'ButtonPushedFcn',@(but,event) buttonpushed(but),'Text','Threshold','UserData',0);
% [center,radius]=imfindcircles(imboxfilt(double(edge((uint16(imgaussfilt(template,3))),'Canny',[0 15/intensityrange(2)],2)),3)>0,[15 100],'Sensitivity',.8,'ObjectPolarity','bright');
% 
% % [center,radius]=imfindcircles(imboxfilt(double(edge(template,'Canny',[0 15/intensityrange(2)],.5)),3)>0,[15 100],'Sensitivity',.8,'ObjectPolarity','bright');
% sens=.85;
% while 1
%     [center,radius]=imfindcircles(imboxfilt(double(edge((uint16(imgaussfilt(template,3))),'Canny',[0 15/intensityrange(2)],2)),3)>0,[15 100],'Sensitivity',sens,'ObjectPolarity','bright');
% if ~isempty(center)
%     break;
% elseif sens==.95
%     break
% else
%     sens=sens+.05;
% end
% end
% if isempty(center)
% %     [center,radius]=imfindcircles(imboxfilt(double(edge(imresize(template,.5),'Canny',[0 15/intensityrange(2)],.5)),1)>0,[10 50],'Sensitivity',.8,'ObjectPolarity','bright');
% [center,radius]=imfindcircles(imboxfilt(double(edge((uint16(imgaussfilt(imresize(template,.5),3))),'Canny',[0 15/intensityrange(2)],2)),3)>0,[15 100],'Sensitivity',.75,'ObjectPolarity','bright');
% 
%     center=center*2;
%     radius=radius*2;
% end
% if isempty(center)
%     [center,radius]=imfindcircles(imboxfilt(double(edge(imresize(template,.5),'Canny',[0 4/intensityrange(2)],.5)),3)>0,[10 50],'Sensitivity',.8,'ObjectPolarity','bright');
%     center=center*2;
%     radius=radius*2;
% end

%%%find pupil

A=(edge(imresize(imgaussfilt(template,4),.5),'Canny',.05,.1));
[center,radius]=imfindcircles(imgaussfilt(double(A),.15)>0,[10 60],"ObjectPolarity",'bright','Sensitivity',.9);
center=center*2;
radius=radius*2;
if isempty(center)
se = strel('disk',5);
bg=imopen((uint16(template)),se);
A=(imboxfilt(double( ((uint16(template))-min(bg(:)))<25),51)>0);
CC=bwconncomp(A);
numPixels = cellfun(@numel,CC.PixelIdxList);
% numPixels(numPixels==max(numPixels))=NaN;
numPixels(numPixels>1e5)=NaN;
[~,pupilcluster]=nanmax(numPixels);
% [xx,yy]=ind2sub(size(A),CC.PixelIdxList{pupilcluster});
STATS=regionprops(CC,{'Centroid','MajorAxisLength','MinorAxisLength','Circularity'});
center=ceil(STATS(pupilcluster).Centroid);
radius=mean([STATS(pupilcluster).MajorAxisLength-101 STATS(pupilcluster).MinorAxisLength-51]);
end
%%%
if ~isempty(center)
    center=center(1,:);
    radius=radius(1);
% position=round([center(1)-radius*3.25,center(2)-radius*2.5,radius*6.5,radius*5]);
rad_fact=max(round(radius/30),1);
position=round([center(1)-20*3.25*rad_fact,center(2)-20*2*rad_fact,20*6.5*rad_fact,20*4*rad_fact]);
position=max(position,1);
try
ax.UserData=uint16(template(floor(position(2):sum(position([2 4]))),floor(position(1):sum(position([1 3])))));
catch
    position=[100 100 100 100];
    ax.UserData=uint16(template(floor(position(2):sum(position([2 4]))),floor(position(1):sum(position([1 3])))));
end
imshow((ax.UserData),'parent',ax,'DisplayRange',intensityrange);
else
%     ax.UserData=uint16(template);
    position=[];
end

slider.Visible='off';button.Visible='off';slider2.Visible='off';
drawnow;
while 1
%     slider.Visible='off';button.Visible='off';
%     ax.UserData=template;
% % imshow(template);
% imshow(template,'parent',ax);
% while 1
% try
% h=drawrectangle(ax);
% position=floor(h.Position);
% ax.UserData=template(floor(position(2):sum(position([2 4]))),floor(position(1):sum(position([1 3]))));
% break;
% end
% end
if ~isempty(position)
y = inputdlg('Window Correct? (y/n)');
if iscell(y) && ~isempty(y)
y=y{1};
end
if ~isempty(y)
switch y
    case 'y'
        break;
    otherwise
%         button.UserData=0;
end
else
    break
end
end
%     slider.Visible='off';button.Visible='off';
    ax.UserData=uint16(template);
% imshow(template);
imshow((ax.UserData),'parent',ax,'DisplayRange',[intensityrange(1) intensityrange(2)]);
drawnow;
while 1
try
h=drawrectangle(ax);
position=floor(h.Position);
ax.UserData=uint16(template(floor(position(2):sum(position([2 4]))),floor(position(1):sum(position([1 3])))));
imshow((ax.UserData),'parent',ax,'DisplayRange',intensityrange);
break;
end
end
end
th=max(ceil(quantile(double(ax.UserData(:)),.01)*1.5),1);
slider.Value=th;
% colorpupil=repmat(ax.UserData,[1,1,3]);
temp=min(ax.UserData,intensityrange(2));
colorpupil=repmat(temp,[1,1,3]);
% colorpupil(repmat(ax.UserData<=th,[1,1,3]))=reshape(repmat([intensityrange(2) intensityrange(2) 0],[sum(ax.UserData(:)<=th),1,1]),[],1);
colorpupil=recolor_pupil(colorpupil,intensityrange,temp,th);
% imshow(double(colorpupil)./intensityrange(2),'parent',ax,'DisplayRange',intensityrange);
imshow(double(colorpupil)./intensityrange(2),'parent',ax,'DisplayRange',[0 1]);
% se=strel('disk',2);
% A=edge(imgaussfilt(double(imopen(ax.UserData<=th,se)),2));
% 
% % [center,radius]=imfindcircles(imresize(imgaussfilt(double(A),2),.5)>0,[10 60],'ObjectPolarity','bright','sensitivity',.8);center=center*2;radius=radius*2;
% % if ~isempty(center) && size(center,1)==1
% %     coords1=repmat(1:size(A,2),size(A,1),1);
% %     coords2=repmat((1:size(A,1))',1,size(A,2));
% %     r2=(coords1-center(1)).^2+(coords2-center(2)).^2;
% %     A(r2<(radius-6).^2)=0;
% % end
% CC=bwconncomp(A);
% if length(CC)>1
% numPixels=cellfun(@numel,CC.PixelIdxList);
% [~,in_max]=max(numPixels);
% [ii,jj]=ind2sub(size(A),CC.PixelIdxList{in_max});
% else
%     [ii,jj]=find(A>0);
% end
% 
% % [ii,jj]=find(A);
% % 
% params = FitEllipse(ii,jj);
% % 
% toexclude=false(length(jj),1);
% for a=1:3
% %     toexclude=(sqrt((jj-params.yc).^2+(ii-params.xc).^2))<min(params.ra,params.rb)-2 | toexclude;
% % xydist=abs((jj-params.yc).^2/params.ra.^2+(ii-params.xc).^2/params.rb.^2-1);
% % toexclude=toexclude | xydist>mean(xydist)+std(xydist);\
% toexclude=toexclude | (sqrt((jj-params.xc).^2+(ii-params.yc).^2))<min(params.ra,params.rb);
% params=FitEllipse(ii(~toexclude),jj(~toexclude));
% end
% slider2.Value=round(rad2deg(params.ang/45)*100);
% 
%         tform = maketform('affine',[1 0 0; slider2.Value/100 1 0; 0 0 1]);
%         colorpupil = imtransform(colorpupil,tform,'bilinear','XData',[1 size(colorpupil,2)],'YData',[1 size(colorpupil,1)]);
% 

% keyboard

while 1
y = inputdlg('Threshold Correct? (y/n)');
if iscell(y) && ~isempty(y)
y=y{1};
end
if ~isempty(y)
switch y
    case 'y'
        break;
    otherwise
        button.UserData=0;
end
else
    break
end

slider.Visible='on';button.Visible='on';slider2.Visible='on';
drawnow;
while button.UserData==0
    drawnow;
end
end
button.Enable='off';
slider.Enable='off';
slider2.Enable='off';
th=slider.Value;
temp=ax.UserData<=th;
        tform = maketform('affine',[1 0 0; slider2.Value/100 1 0; 0 0 1]);
        temp = imtransform(temp,tform,'bilinear','XData',[1 size(temp,2)],'YData',[1 size(temp,1)]);
se=strel('disk',1);
A=edge(imgaussfilt(double(imopen(temp,se)),2));

% [center,radius]=imfindcircles(imresize(imgaussfilt(double(A),2),.5)>0,[10 60],'ObjectPolarity','bright','sensitivity',.8);center=center*2;radius=radius*2;
% if ~isempty(center) && size(center,1)==1
%     coords1=repmat(1:size(A,2),size(A,1),1);
%     coords2=repmat((1:size(A,1))',1,size(A,2));
%     r2=(coords1-center(1)).^2+(coords2-center(2)).^2;
%     A(r2<(radius-6).^2)=0;
% end
CC=bwconncomp(A);
if length(CC)>1
numPixels=cellfun(@numel,CC.PixelIdxList);
[~,in_max]=max(numPixels);
[ii,jj]=ind2sub(size(A),CC.PixelIdxList{in_max});
else
    [ii,jj]=find(A>0);
end

% [ii,jj]=find(A);
% 
try
params = FitEllipse(ii,jj);
% 
toexclude=false(length(jj),1);
for a=1:3
%     toexclude=(sqrt((jj-params.yc).^2+(ii-params.xc).^2))<min(params.ra,params.rb)-2 | toexclude;
% xydist=abs((jj-params.yc).^2/params.ra.^2+(ii-params.xc).^2/params.rb.^2-1);
% toexclude=toexclude | xydist>mean(xydist)+std(xydist);\
toexclude=toexclude | (sqrt((jj-params.xc).^2+(ii-params.yc).^2))<min(params.ra,params.rb);
params=FitEllipse(ii(~toexclude),jj(~toexclude));
end

% % params.xc^2/params.ra^2+params.yc^2/params.rb^2=1;
% outlie=false(length(ii),1);
% for a=1:5
% outlie=outlie | abs((jj-params.xc).^2/params.ra.^2+(ii-params.yc).^2/params.rb.^2-1)>(.1-(a-1)*.01);
% try
% params = FitEllipse(ii(~outlie),jj(~outlie));
% end
% end
position_p=[params.xc-params.ra params.yc-params.rb params.ra*2 params.rb*2];
catch
    position_p=[100 100 100 100];
end
h3=rectangle('Position',position_p,'curvature',1,'LineStyle','-','EdgeColor','r','LineWidth',3,'Parent',ax);
while 1
    y = inputdlg('Pupil Correct? (y/n)');
if iscell(y) && ~isempty(y)
y=y{1};
end
if ~isempty(y)
switch y
    case 'y'
        break;
    otherwise
%         button.UserData=0;
delete(h3);
end
else
    break
end
while 1
try
h2=drawrectangle(ax);
position_p=floor(h2.Position);
delete(h2);
h3=rectangle('Position',position_p,'curvature',1,'LineStyle','-','EdgeColor','r','LineWidth',3,'Parent',ax);
break;
end
end


end

% temp=ax.UserData;
temp=min(ax.UserData,intensityrange(2));
colorpupil=repmat(temp,[1,1,3]);
% colorpupil(repmat(ax.UserData<=th,[1,1,3]))=reshape(repmat([intensityrange(2) intensityrange(2) 0],[sum(ax.UserData(:)<=th),1,1]),[],1);
colorpupil=recolor_pupil(colorpupil,intensityrange,temp,th);
temp=colorpupil;

        tform = maketform('affine',[1 0 0; slider2.Value/100 1 0; 0 0 1]);
        temp = imtransform(temp,tform,'bilinear','XData',[1 size(temp,2)],'YData',[1 size(temp,1)]);
A=edge(imgaussfilt(double(ax.UserData),2),'Canny',.06,1);
A(floor(position_p(2))-2:ceil(position_p(4)+position_p(2))+2,floor(position_p(1))-2:ceil(position_p(3)+position_p(1))+2)=0;
BW=1-(imboxfilt(double(A),9)>0);
CC=bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);[~,in]=max(numPixels);for a=setdiff(1:length(numPixels),in);BW(CC.PixelIdxList{a})=0;end
se=strel('disk',40);
BW=(imclose(BW,se))>0;
temp(~BW)=intensityrange(2);
temp(1,:)=intensityrange(2);
temp(end,:)=intensityrange(2);
temp(:,1)=intensityrange(2);
temp(:,end)=intensityrange(2);
% imshow(min(temp*3,intensityrange(2)),'parent',ax,'DisplayRange',intensityrange);
imshow(double(temp)./intensityrange(2),'parent',ax,'DisplayRange',[0 1]);

while 1
%     slider.Visible='off';button.Visible='off';
%     ax.UserData=template;
% % imshow(template);
% imshow(template,'parent',ax);
% while 1
% try
% h=drawrectangle(ax);
% position=floor(h.Position);
% ax.UserData=template(floor(position(2):sum(position([2 4]))),floor(position(1):sum(position([1 3]))));
% break;
% end
% end
if ~isempty(position)
y = inputdlg('Background Correct? (y/n)');
if iscell(y) && ~isempty(y)
y=y{1};
end
if ~isempty(y)
switch y
    case 'y'
        break;
    otherwise
%         button.UserData=0;
end
else
    break
end
end
%     slider.Visible='off';button.Visible='off';
%     ax.UserData=uint16(template);
% temp=ax.UserData;
temp=min(ax.UserData,intensityrange(2));
colorpupil=repmat(temp,[1,1,3]);
% colorpupil(repmat(ax.UserData<=th,[1,1,3]))=reshape(repmat([intensityrange(2) intensityrange(2) 0],[sum(ax.UserData(:)<=th),1,1]),[],1);
colorpupil=recolor_pupil(colorpupil,intensityrange,temp,th);
temp=colorpupil;
        tform = maketform('affine',[1 0 0; slider2.Value/100 1 0; 0 0 1]);
        temp = imtransform(temp,tform,'bilinear','XData',[1 size(temp,2)],'YData',[1 size(temp,1)]);
% imshow(min(temp*3,intensityrange(2)),'parent',ax,'DisplayRange',[0 intensityrange(2)]);
imshow(double(temp)./intensityrange(2),'parent',ax,'DisplayRange',[0 1]);


% imshow(template);
% imshow((ax.UserData),'parent',ax,'DisplayRange',[0 intensityrange(2)]);
drawnow;
while 1
try
h_p=drawpolygon(ax);
position_pol=floor(h_p.Position);
delete(h_p)
non_bg=inpolygon(repmat(1:size(ax.UserData,2),size(ax.UserData,1),1),repmat((1:size(ax.UserData,1))',1,size(ax.UserData,2)),position_pol(:,1),position_pol(:,2));
bg=imgaussfilt(double(~non_bg),1,'Padding',NaN);
bg(isnan(bg))=1;
bg=bg>0;
temp(bg)=intensityrange(2);
% imshow(min(temp*3,intensityrange(2)),'parent',ax,'DisplayRange',[0 intensityrange(2)]);
imshow(double(temp)./intensityrange(2),'parent',ax,'DisplayRange',[0 1]);
BW=~bg;
break;
end
end
end



pupil_aspectratio=position_p(3)./position_p(4);
pupil_radius=max(position_p(3),position_p(4));

mv=slider.Value;
skew_val=slider2.Value;

input=struct;
input.vid_file=vid_file;
input.eye_outline=BW;
input.skew=skew_val;
input.eye_window=[position(1) sum(position([1,3])) position(2) sum(position([2,4]))];
if exist('vidObj')
input.num_frames=floor(vidObj.Duration*60*vidObj.FrameRate);
else
    input.num_frames=NaN;
end
input.mv=mv;
input.pupil_aspectratio=pupil_aspectratio;
input.pupil_radius=pupil_radius;
input.params=params;
input.se=strel('disk',2);
input.coords_y=repmat((1:size(temp,1))',1,size(temp,2));
input.coords_x=repmat((1:size(temp,2)),size(temp,1),1);
close(f);
% keyboard

end

function [event_mat,rot_val,timestamps,events]=read_treadmill_timestamps(TS_file)
fid=fopen(TS_file,'r');
[data]=fscanf(fid,'%u%u%u%u%u%u%u%u,%d,%ld');
data=reshape(data,10,[])';
rot_val=data(:,9);
maxval=max(rot_val);
transitions=cumsum([0;(diff(rot_val)<-100)].*(maxval+1));
rot_val=rot_val+transitions;
timestamps=data(:,10);
event_mat=logical(data(:,1:8));
events=struct();
struct_fields={'pupil','twophoton','lick_water','lick_no_water','sound','water_no_lick','shock','movement'};
for struct_rep=1:length(struct_fields)
setfield(events,struct_fields{struct_rep},timestamps(event_mat(:,struct_rep)));
end
end

function updateColorMap(sld2,sld,ax,intensityrange)
%     caxis(ax,[0 sld.value]);
th=sld2.Value;
% imshow(ax.UserData>sld.Value,'parent',ax,'DisplayRange',[0 1]);
temp=min(ax.UserData,intensityrange(2));
% temp=imrotate(temp,sld.Value,'bilinear','crop');
% temp=imrotate(temp,sld.Value,'bilinear','crop');
        tform = maketform('affine',[1 0 0; sld.Value/100 1 0; 0 0 1]);
        temp = imtransform(temp,tform,'bilinear','XData',[1 size(temp,2)],'YData',[1 size(temp,1)]);

% temp(temp==0)=intensityrange(2);
colorpupil=repmat(temp,[1,1,3]);
% colorpupil(colorpupil<=th*3)=reshape(repmat([intensityrange(2) intensityrange(2) 0],[sum(temp(:)<=th*3),1,1]),[],1);
colorpupil=recolor_pupil(colorpupil,intensityrange,temp,th);
% imshow(colorpupil,'parent',ax,'DisplayRange',intensityrange);
imshow(double(colorpupil)./intensityrange(2),'parent',ax,'DisplayRange',[0 1]);

end

function sld_rotate(sld,sld2,ax,intensityrange)
th=sld2.Value;
% imshow(ax.UserData>sld.Value,'parent',ax,'DisplayRange',[0 1]);
temp=min(ax.UserData,intensityrange(2));
% temp=imrotate(temp,sld.Value,'bilinear','crop');
        tform = maketform('affine',[1 0 0; sld.Value/100 1 0; 0 0 1]);
        temp = imtransform(temp,tform,'bilinear','XData',[1 size(temp,2)],'YData',[1 size(temp,1)]);

% temp(temp==0)=intensityrange(2);
colorpupil=repmat(temp,[1,1,3]);
% colorpupil(colorpupil<=th*3)=reshape(repmat([intensityrange(2) intensityrange(2) 0],[sum(temp(:)<=th*3),1,1]),[],1);

colorpupil=recolor_pupil(colorpupil,intensityrange,temp,th);
% imshow(colorpupil,'parent',ax,'DisplayRange',intensityrange);
imshow(double(colorpupil)./intensityrange(2),'parent',ax,'DisplayRange',[0 1]);
end

function buttonpushed(but)
but.UserData=1;
end
%1=pupil camera
%2=microscope camera
%3=lick, water available
%4=lick, water not available
%5=sound
%6=unprompted water
%7=shock
%8=run
%9=timestamp

function out=recolor_pupil(colorpupil,intensityrange,temp,th)
colorpupil(colorpupil<=th)=...
    min(colorpupil(colorpupil<=th).*uint16(reshape(repmat([3 3 1],sum(temp(:)<=th),1),[],1)),intensityrange(2));
zero_fix=sum(sum(colorpupil(:,:,1)==0));
colorpupil(colorpupil==0)=[intensityrange(2)*ones(zero_fix,1,'uint16');intensityrange(2)*ones(zero_fix,1,'uint16');zeros(zero_fix,1,'uint16')];
out=colorpupil;
end