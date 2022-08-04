function data_out=RR_from_thermal(data_out)
if isfield(data_out,'thermal_data')
    for a=1:length(data_out.thermal_data)
        peaks=data_out.thermal_data(a).inhale_stats.starts_idx(data_out.thermal_data(a).inhale_stats.length>30000);
 RR=1e6*60*9./[double(data_out.thermal_data(a).ts(peaks(10:end)))-double(data_out.thermal_data(a).ts(peaks(1:end-9)))];
%             data_out.respiration_peaks=peaks(:);
            peaktimes=peaks(5:end-5);
            RR_full=interp1(peaktimes,RR,1:length(data_out.thermal_data(a).ts),'linear');
            RR_full(isnan(RR_full))=interp1(find(~isnan(RR_full)),RR_full(~isnan(RR_full)),find(isnan(RR_full)),'nearest','extrap');
            data_out.thermal_data(a).RR=RR_full(:);

            binary_mat=nan(1,length(data_out.thermal_data(a).RR));
            % peaktimes(tofix)=NaN;
            mid_point_pt=min(max(round((peaktimes(2:end)+peaktimes(1:end-1))/2),1),length(binary_mat));
            binary_mat(mid_point_pt)=diff(data_out.thermal_data(a).ts(peaktimes)/1e3);
            binary_mat=binary_mat./nanmean(binary_mat);
            RRV=sqrt(movmean(binary_mat.^2,1*100,'omitnan'));

            RRV_interp=medfilt1(RRV,3*100,'omitnan');
            data_out.thermal_data(a).peak_vector=binary_mat;
            RRV(isnan(RRV))=interp1(find(~isnan(RRV)),RRV_interp(~isnan(RRV)),find(isnan(RRV)),'linear');
            RRV(isnan(RRV))=interp1(find(~isnan(RRV)),RRV(~isnan(RRV)),find(isnan(RRV)),'nearest','extrap');
            RRV=movmean(RRV,30);
            data_out.thermal_data(a).RRV=RRV(:);
    end
end