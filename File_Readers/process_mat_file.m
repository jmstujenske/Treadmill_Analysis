function [data_out]=process_mat_file(data_out);
     temp=load(data_out.mat_path);
     if ~isfield(temp,'Laser')
         temp.Laser.time=[];
         temp.Laser.duration=[];
     end
    matdata=temp;
    fulltime=data_out.timestamps(end)-data_out.timestamps(1);
    
    real_mattimes=matdata.Events.time-matdata.Events.time(1)+matdata.options.FirstTone_delay*1e3;
    todelete=real_mattimes>fulltime/1e3;
    S=fieldnames(matdata.Events);
    if any(todelete)
    for nS=1:length(S)
        if ~all(todelete)
        matdata.Events.(S{nS})=matdata.Events.(S{nS})(~todelete);
        else
            matdata.Events.(S{nS})=[];
        end
    end
    end
    real_mattimes=matdata.Sounds.time-matdata.Sounds.time(1)+matdata.options.FirstTone_delay*1e3;
    todelete=real_mattimes>fulltime;
    S=fieldnames(matdata.Sounds);
    if any(todelete)
    for nS=1:length(S)
        if ~all(todelete)
        matdata.Sounds.(S{nS})=matdata.Sounds.(S{nS})(~todelete);
        else
            matdata.Sounds.(S{nS})=[];
        end
    end
    end
    n_event=length(matdata.Events.duration);
    events.types=zeros(1,n_event);
        events.types(matdata.Events.duration>1000)=1;
        if isfield(matdata.Events,'water')
            confield='water';
        elseif isfield(matdata.Events,'shock')
            confield='shock';
        end
%     events.types(matdata.Events.shock>0)=2;
    events.types(matdata.Events.laser>0)=3;
    n_ts=length(matdata.Events.(confield));
    data_out.sound_times=data_out.event_times(events.types(1:n_ts)==1);
    data_out.sound_frequency=matdata.Sounds.frequency;
    data_out.laser_times=data_out.event_times(events.types(1:n_ts)==3);
shockorlaser=cumsum(matdata.Events.(confield)>0 | matdata.Events.laser>0);
lasertones=(find(matdata.Events.laser>0)-shockorlaser(matdata.Events.laser>0)+1);
data_out.sound_withlaser=lasertones;