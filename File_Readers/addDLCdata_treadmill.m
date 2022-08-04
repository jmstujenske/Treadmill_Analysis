function data_out=addDLCdata_treadmill(data_out,video_dir)
            if ~isempty(data_out.bin_path)
%             data_out.resp_path=csv_file;
            if isfield(data_out,'resp_path') && ~isempty(data_out.resp_path)
            csv_file=data_out.resp_path;
            else
                                if ~isstruct(data_out.bin_path)
            ins=regexp(data_out.bin_path,'\');
            animal_name=data_out.bin_path(ins(end)+1:end-10);
                else
                    animal_name=data_out.bin_path(1).name(1:end-10);
                end
            DLC_output=dir([video_dir,animal_name,'*.csv']);
            csv_file=dir2file(DLC_output);
            if strcmp(csv_file,'\')
                csv_file=[];
            end
            data_out.animal_name=animal_name;
            data_out.resp_path=csv_file;
            end
            
        if ~isempty(csv_file)
            
            DLC_data=importdata(csv_file);
            npoints=(size(DLC_data.data,2)-1)/3;
            coords=[];
            for rep=1:npoints
                coords=cat(3,coords,DLC_data.data(:,(2:4)+(rep-1)*3));
            end
            coords2=coords;
            for rep=1:size(coords,3)
                coords2(coords2(:,3,rep)<.9,1:2,rep)=NaN;
                cnt=hist3(coords2(:,1:2,rep),{1:960,1:600});
                cnt=cnt';
                vals_in=find(imgaussfilt(double(cnt./sum(cnt(:))>.0005),10)>0);
                inds2=sub2ind([600 960],nanmax(nanmin(round(coords(:,2,rep)),600),1),nanmax(nanmin(round(coords(:,1,rep)),960),1));
                coords2(~ismember(inds2,vals_in) & coords2(:,3,rep)<.99,1:2,rep)=NaN;
                indices=1:length(coords2);
%                 errors=false(1,sum(~isnan(coords2(:,1,rep)))-1);
                while 1
                    notnan=find(~isnan(coords2(:,1,rep)));
                    errors=abs(diff(coords2(notnan,1,rep))./diff(notnan))>20;
                    coords2(notnan(find(errors)+1),1:2,rep)=NaN;
                    if ~any(errors)
                        break
                    end
                end
                coords2(coords2(:,3,rep)<.95,1:2,rep)=NaN;
                tofix=isnan(coords2(:,1,rep));
                for xy_rep=1:2
                    coords2(:,xy_rep,rep)=interp1(indices(~tofix),coords2(~tofix,xy_rep,rep),indices,'linear','extrap');
                end
            end
            coords2=medfilt3(coords2,[5 1 1]);
%             coords2=coords2-nanmean(coords2(:,:,8:10),3);
            coords_pca=nan(size(coords2,1),9);
            for rep=1:size(coords,3)
                [coeff,score]=pca([coords2(:,1:2,rep)]);
                coords_pca(:,rep)=score(:,1);
            end
            [U S V]=svdecon(coords_pca(:,[1:7]));
            
            d1=zscore(U(:,1));
%             d1=nanmean(zscore(coords_pca(:,[2:3 5:7])),2);
            
            signal=movmax(movmin(LowFilt_Order(d1,30,30,1),30),30);
            sniff_start=conv(diff(signal)<-.035,ones(1,3),'full');
            sniff_start=find(sniff_start(1:end-1));
            sniff_stop=find(diff(signal)>.035);
            %resp=BandFilt_Order(d-movmean(d,30),30,30,1,10);
            vid_times=double(data_out.video_frame_times(1:length(d1)));
            resp=d1;
            upscale_fact=1;
            data_out.sniff_up=sniff_start;
            data_out.sniff_down=sniff_stop;
            data_out.nose_coordinates=coords2;
            data_out.nose_coordinates_pca=coords_pca;
            data_out.resp_signal=resp;
            [m,~,m_pow]=BandFilt_Order(d1,30,30,2,6);
            [m2,~,m2_pow]=BandFilt_Order(d1,30,30,5,9);
            m=m.*(rms(m2)./rms(m));
%             m_pow=m_pow.*(rms(m2)./rms(m));
m_pow=movmin(movmax(m,9),9);
m2_pow=movmin(movmax(m2,7),7);
            low_segs=m_pow>=m2_pow;
            high_segs=m2_pow>m_pow;
            [peaks1,troughs1,p_amp,t_amp]=findpeaks(m);
            [peaks2,troughs2,p_amp,t_amp]=findpeaks(m2);
            concat=zeros(size(m));
            concat(peaks1)=1;
            concat(peaks2)=2;
            concat(intersect(peaks1,peaks2))=3;
            concat(concat==1 & high_segs)=0;
            concat(concat==2 & low_segs)=0;
            data_out.resp_filtered_low=m;
            data_out.resp_filtered_high=m2;
            peaks=find(concat>0);
            peaks(find(diff(peaks)<=2)+1)=[];
            
%             th=quantile(p_amp,.01);
%             if peaks1(1)>troughs1(1)
%                 toexclude=union(find(abs(t_amp)<th),find(p_amp<th));
%             else
%                 toexclude=union(find(abs(t_amp)<th)-1,find(p_amp<th));
%                 toexclude(toexclude==0)=[];
%             end
%             toexclude=intersect(toexclude,[find(diff(peaks1)<=5) find(diff(peaks1)<=5)+1]);
%             peaks1(toexclude)=[];
            RR=1e6*60*9./[double(vid_times(ceil(peaks(10:end)/upscale_fact)))-double(vid_times(ceil(peaks(1:end-9)/upscale_fact)))];
            data_out.respiration_peaks=peaks(:);
            peaktimes=peaks(5:end-5);
            RR_full=interp1(peaktimes,RR,1:length(resp),'linear');
            RR_full(isnan(RR_full))=interp1(find(~isnan(RR_full)),RR_full(~isnan(RR_full)),find(isnan(RR_full)),'nearest','extrap');
            data_out.RR=RR_full(:);
            [yo, fo, to]=mtcsg(d1-movmean(d1,30),30,30,30,29,1);
            [~,in]=max((yo(fo<=10,:).*(fo(fo<=10).^2))./sum(yo(fo<=10,:).*(fo(fo<=10).^2)));
            data_out.RR2=fo(in).*60;
            data_out.RR2=[data_out.RR2(1)*ones(15,1);data_out.RR2;data_out.RR2(end)*ones(15,1)];
            data_out.RR2=data_out.RR2(:);
            binary_mat=nan(1,length(data_out.RR));
            % peaktimes(tofix)=NaN;
            mid_point_pt=min(max(round((peaktimes(2:end)+peaktimes(1:end-1))/2),1),length(binary_mat));
            binary_mat(mid_point_pt)=diff(vid_times(peaktimes)/1e3);
            binary_mat=binary_mat./nanmean(binary_mat);
            RRV=sqrt(movmean(binary_mat.^2,10*30,'omitnan'));

            RRV_interp=medfilt1(RRV,30*30,'omitnan');
            data_out.peak_vector=binary_mat;
            RRV(isnan(RRV))=interp1(find(~isnan(RRV)),RRV_interp(~isnan(RRV)),find(isnan(RRV)),'linear');
            RRV(isnan(RRV))=interp1(find(~isnan(RRV)),RRV(~isnan(RRV)),find(isnan(RRV)),'nearest','extrap');
            RRV=movmean(RRV,9);
            data_out.RRV=RRV(:);
        end
            else
                data_out.resp_path=[];
            end