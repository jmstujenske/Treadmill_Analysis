function data_out=generate_struct_treadmill(animal_name,time_dir,video_dir,mat_dir,thermal_dir)
if nargin<3 || isempty(video_dir)
    video_dir=time_dir;
end
if nargin<4 || isempty(mat_dir)
    mat_dir=time_dir;
end
if nargin<5 || isempty(thermal_dir)
    thermal_dir=video_dir;
end
timestamp_file=['*',animal_name,'*.txt'];
        bin_file=['*',animal_name,'*_pupil.bin'];
        matlab_file=['*',animal_name,'*.mat'];
        timestamp_path=dir([time_dir,timestamp_file]);
            if length(timestamp_path)==1
                timestamp_path=[timestamp_path.folder,'\',timestamp_path.name];
            end
            bin_path=dir([video_dir,bin_file]);
            if length(bin_path)>1
                bin_file=['*',animal_name,'*_pupil.bin'];
                bin_path=dir([video_dir,bin_file]);
            end
            bin_path=[bin_path.folder,'\',bin_path.name];
            if isempty(mat_dir)
                mat_path=[];
            else
            mat_path=dir([mat_dir,matlab_file]);
            end
            if ~isempty(mat_path)
%                 matlab_file=['*',animal_name,'*.mat'];
%                 mat_path=dir([mat_dir,matlab_file]);
                mat_path=[mat_path(end).folder,'\',mat_path(end).name];
%             else
%                mat_path=[mat_path(end).folder,'\',mat_path(end).name];
% 
%             end
            end
            clear data_out
            try
            data_out=extractdata_treadmill_timestamps_HR(timestamp_path);
            end
            data_out.bin_path=bin_path;
            data_out.mat_path=mat_path;
            data_out.timestamp_path=timestamp_path;
            try
            data_out=extract_pupildata(data_out);
            end
            try
            data_out=process_mat_file(data_out);
            end
            data_out.animal_name=animal_name;
            try
            data_out=treadmill_processmovementandlicks(data_out);
            end
            data_out=addDLCdata_treadmill(data_out,video_dir);
%             if ~isempty(data_out.HR)
            try
            data_out=HR_fix2(data_out);
            end
            try
                thermal_files=dir([thermal_dir,'*',animal_name,'*_thermal*']);
                thermal_data=process_thermal_treadmill(thermal_files,data_out.thermal_cam_ontime);
                data_out.thermal_data=thermal_data;
                data_out=RR_from_thermal(data_out);
            end