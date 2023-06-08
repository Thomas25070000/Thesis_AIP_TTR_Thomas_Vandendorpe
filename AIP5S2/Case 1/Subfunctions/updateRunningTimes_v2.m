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
        block_sections = timetable.blocksection(trains(tt).rows);
        if ismember(trains(tt).type,{'IC','THA'})
            if trains(tt).direction == 2
                RT = runningtimes.IC02.disrupted(block_sections);
            end
            if trains(tt).direction == 3
                RT = runningtimes.IC03.disrupted(block_sections);
            end
            if trains(tt).direction == 12
                RT = runningtimes.IC12.disrupted(block_sections);
            end
            if trains(tt).direction == 13
                RT = runningtimes.IC13.disrupted(block_sections);
            end
        else
            if trains(tt).direction == 2
                RT = runningtimes.L02.disrupted(block_sections);
            end
            if trains(tt).direction == 3
                RT = runningtimes.L03.disrupted(block_sections);
            end
            if trains(tt).direction == 12
                RT = runningtimes.L12.disrupted(block_sections);
            end
            if trains(tt).direction == 13
                RT = runningtimes.L13.disrupted(block_sections);
            end
        end
        timetable.running(trains(tt).rows) = RT;
    end
end


    
    
    
    
end