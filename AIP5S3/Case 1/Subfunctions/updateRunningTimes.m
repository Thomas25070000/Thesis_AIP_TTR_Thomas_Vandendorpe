function timetable = updateRunningTimes(timetable, blocksections, settings)

v_S1 = settings.infrastructure.switch.maxvS1 / 3.6;
v_S2 = settings.infrastructure.switch.maxvS2 / 3.6;
v_closed = settings.disruption.maxspeed / 3.6;

v_S1 = min(v_S1,v_closed);
v_S2 = min(v_S2,v_closed);

for ee = 1:size(timetable,1)
    bs = timetable.blocksection(ee);
    l_block = blocksections.length(bs);
    
    try
        if timetable.train_id(ee) ~= timetable.train_id(ee-1)
            first_ev_train = ee;
        end
    catch
        first_ev_train = 1;
    end
    
    if timetable.arrival(first_ev_train) < settings.disruption.duration*3600
        time = timetable.running(ee);
        switch blocksections.type{bs}
            case 'S1'
                if blocksections.closed(bs)
                    % Speed has to be reduced, so time increased.
                    time = round(l_block / v_S1);
                end
            case 'D'
                time = round(l_block / v_closed);
            case 'S2'
                if blocksections.closed(bs)
                    % Speed has to be reduced, so time increased.
                    time = round(l_block / v_S2);
                end
        end
        timetable.running(ee) = time;
    end
end