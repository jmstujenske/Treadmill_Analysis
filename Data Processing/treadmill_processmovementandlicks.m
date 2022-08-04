function data_struct=treadmill_processmovementandlicks(data_struct)
allpos_change=diff(double(data_struct.position));
allpos_change(allpos_change<-100)=1;
allpos_change(allpos_change>100)=-1;
allpos_change(allpos_change<-3)=0;
allpos_change(allpos_change>3)=0;
newpos=cumsum([0;allpos_change]);
newpos=movmedian(movmin(movmax(movmedian(newpos,21),1001),1001),3001);
% nomotion=movmax(movmin([diff(newpos)==0;true],51),51);
% newpos=movmean(newpos,21);
newpos=LowFilt_Order(newpos,1000,10,5);
newpos(end-25:end)=newpos(end-25);
velocity=double(newpos(21:end)-newpos(1:end-20))./double(data_struct.timestamps(21:end)-data_struct.timestamps(1:end-20));

velocity=[velocity(1)*ones(10,1);velocity;velocity(end)*ones(10,1)];
% velocity(nomotion)=0;
data_struct.velocity=velocity;
data_struct.position_linear=newpos;
temp=double(data_struct.licks_times);

       boutends=find(diff(temp)>250e3);
        boutstarts=[1;boutends+1];
        boutends=[boutends;length(temp)];
realbouts=boutends-boutstarts>2;
        data_struct.lickbouts_times=[temp(boutstarts(realbouts)) temp(boutends(realbouts))];
        data_struct.licks_times_within_bout=[];
        data_struct.boutlicking_logical=false(1,length(data_struct.timestamps));

        for b=1:sum(realbouts)
            thisbout_data=data_struct.licks_times(temp>=data_struct.lickbouts_times(b,1) & temp<=data_struct.lickbouts_times(b,2));
            data_struct.licks_times_within_bout=[data_struct.licks_times_within_bout;thisbout_data];
            times=data_struct.timestamps>=data_struct.lickbouts_times(b,1) & data_struct.timestamps<=data_struct.lickbouts_times(b,2);
            data_struct.boutlicking_logical(times)=true;
        end
