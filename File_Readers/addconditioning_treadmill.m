function data_struct=addconditioning_treadmill(data_struct)
%days=first day of shock

for animal_rep=1:size(data_struct,2)
    first_reward=true;
    first_shock=true;
    for days_rep=1:size(data_struct,1)
        try
        matdata=load(data_struct(days_rep,animal_rep).mat_path);
        if isfield(matdata.Events,'water')
            data_struct(days_rep,animal_rep).Conditioning_type='Reward';
            if first_reward
                reward_tone=matdata.Sounds.frequency(find(matdata.Sounds.water,1,'first'));
                first_reward=false;
            end
            data_struct(days_rep,animal_rep).reward_freq=reward_tone;
            data_struct(days_rep,animal_rep).sounds_CSplus=data_struct(days_rep,animal_rep).sound_frequency==reward_tone;
            data_struct(days_rep,animal_rep).sound_withreward=matdata.Sounds.water;
        elseif isfield(matdata.Events,'shock')
            data_struct(days_rep,animal_rep).Conditioning_type='Fear';
            
            if first_shock
                for a=1:size(data_struct,1)
                    if ~isempty(data_struct(a,animal_rep).shock_times)
                        [~,in]=min(abs((double(data_struct(a,animal_rep).sound_times)-double(data_struct(a,animal_rep).shock_times(1)))/1e6+29));
                        shock_tone=data_struct(a,animal_rep).sound_frequency(in);
                        break;
                    end
                end
                first_shock=false;
            end
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
    end
end