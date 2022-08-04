function data_struct=addshockdata_treadmill(data_struct,day)
%days=first day of shock
if nargin<2 || isempty(day)
    day=find(~cellfun(@isempty,{data_struct(:,1).shock_times}),1,'first');
end

for animal_rep=1:size(data_struct,2)
    [~,in]=min(abs((double(data_struct(day,animal_rep).sound_times)-double(data_struct(day,animal_rep).shock_times(1)))/1e6+29));
    shock_tone=data_struct(day,animal_rep).sound_frequency(in);
    for days_rep=day:size(data_struct,1)
        data_struct(days_rep,animal_rep).shock_freq=shock_tone;
        data_struct(days_rep,animal_rep).sounds_CSplus=data_struct(days_rep,animal_rep).sound_frequency==shock_tone;
        N_shock=length(data_struct(days_rep,animal_rep).shock_times);
        data_struct(days_rep,animal_rep).sound_withshock=[];
        if N_shock>0
            for shock_rep=1:N_shock
                [~,in]=min(abs((double(data_struct(days_rep,animal_rep).sound_times)-double(data_struct(days_rep,animal_rep).shock_times(shock_rep)))/1e6+29));
                data_struct(days_rep,animal_rep).sound_withshock=[data_struct(days_rep,animal_rep).sound_withshock in];
            end
        end
    end
end