function [center,radius,area]=findpupil_inframe(frame,input,previouscenter,previousradius)
if nargin<3 || isempty(previouscenter)
    previouscenter=[NaN NaN];
end
if nargin<4 || isempty(previousradius)
    previousradius=[NaN NaN];
end
if any(~isnan(previousradius))
radius2=((input.coords_y-previouscenter(2))./previousradius(2)).^2+((input.coords_x-previouscenter(1))./previousradius(1)).^2;
mask=(radius2<1.5 & radius2>.5);
else
    mask=NaN;
end
frame=frame(input.eye_window(3):input.eye_window(4),input.eye_window(1):input.eye_window(2),1);


% if max(frame(:))>1
% %     frame=imgaussfilt(imopen(single(frame<=(input.mv*255+1)),input.se),1,'Padding','symmetric');
% frame=imfilter(imopen(single(frame<=(input.mv*255+1)),input.se), h, 'symmetric', 'conv', 'same');
% else
% % frame=imgaussfilt(imopen(single(frame<=(input.mv+1/255)),input.se),1,'Padding','symmetric');
% frame=imfilter(imopen(single(frame<=(input.mv+1/255)),input.se), h, 'symmetric', 'conv', 'same');
% end
binary_im=false(size(frame));
% filt=imclose(imopen(frame<=(input.mv+1) &input.eye_outline,strel('disk',1)),strel('disk',5));
filt=frame<=(input.mv+1) &input.eye_outline;

% CC=bwconncomp(frame<=(input.mv+1) & input.eye_outline);
CC=bwconncomp(filt);

if CC.NumObjects>0
numPixels=cellfun(@numel,CC.PixelIdxList);
[~,in_max]=max(numPixels);
% [ii,jj]=ind2sub(size(BW1_out),CC.PixelIdxList{in_max});
binary_im(CC.PixelIdxList{in_max})=true;
end
% binary_im=binary_im ;
% h = images.internal.createGaussianKernel([1 1], [7 7]);
% frame=imfilter(imopen(single(binary_im),input.se), h, 'symmetric', 'conv', 'same')>.1;
frame=imfilter(imopen(single(binary_im),input.se), ones(15,15)/15^2, 'symmetric', 'conv', 'same')>.5;
% frame=imfilter(single(frame>=.6),h,'symmetric','conv','same');
% frame=imfilter(single(binary_im), h, 'symmetric', 'conv', 'same')>.6;
try
[center,radius,area]=findpupil(frame,mask);
catch
    center=[NaN NaN];
    radius=[NaN NaN];
    area=NaN;
end
end

function [center,radius,area]=findpupil(G_in,mask)

BW1_out=single(edge(G_in,'Sobel'));
if nargin>2 && ~isempty(mask) && numel(mask)>1
    BW1_out=BW1_out & mask;
end
% BW1_out(~bg)=0;
% CC=bwconncomp(BW1_out);
% if length(CC)>1
% numPixels=cellfun(@numel,CC.PixelIdxList);
% [~,in_max]=max(numPixels);
% [ii,jj]=ind2sub(size(BW1_out),CC.PixelIdxList{in_max});
% else
    [ii,jj]=find(BW1_out>0);
% end
if ~isempty(ii)
try
params=FitEllipse_brief(ii,jj);
toexclude=false(length(jj),1);
n_exclude_steps=3;
for a=1:n_exclude_steps
% toexclude=toexclude | (sqrt((jj-params.xc).^2+(ii-params.yc).^2))<min(params.ra,params.rb)-(3-a);

toexclude=toexclude | sqrt((jj-params.xc).^2+(ii-params.yc).^2./(params.rb^2).*(params.ra^2))<(params.ra-(n_exclude_steps-a));

    params=FitEllipse_brief(ii(~toexclude),jj(~toexclude));
end
center=[params.xc params.yc];
radius=[params.ra params.rb];
area=pi*params.ra*params.rb;
catch
center=[NaN NaN];
radius=[NaN NaN];
area=NaN;
end
else
    center=[NaN NaN];
radius=[NaN NaN];
area=NaN;
end
end