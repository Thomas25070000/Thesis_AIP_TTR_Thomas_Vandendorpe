function timetable = updateRunningTimes_v2(timetable, runningtimes, settings)

train_ids = unique(timetable.train_id);
trains = [];
for tt = 1:length(train_ids)
    rows = find(timetable.train_id == train_ids(tt));
    trains(tt).rows = rows;
    trains(tt).direction = timetable.direction(rows(1));
    trains(tt).type = timetable.train_type(rows(1));
    trains(tt).release = timetable.arrival(rows(1));
end

for tt = 1:length(train_ids)
    
    if trains(tt).release < settings.disruption.duration * 3600
        if ismember(trains(tt).type,{'IC','THA'})
            if trains(tt).direction == 0
                RT = runningtimes.IC0.disrupted(end:-1:1);
            else
                RT = runningtimes.IC1.disrupted;
            end
        else
            if trains(tt).direction == 0
                RT = runningtimes.L0.disrupted(end:-1:1);
            else
                RT = runningtimes.L1.disrupted;
            end
        end
        timetable.running(trains(tt).rows) = RT;
    end
end


    
    
    
    
end