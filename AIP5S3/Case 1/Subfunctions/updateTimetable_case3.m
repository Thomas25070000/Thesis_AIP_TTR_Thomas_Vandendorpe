function new_timetable = updateTimetable_case3(timetable, settings, traininfo)

Ntrains = size(traininfo,2);
Ntracks = size(settings.tracks,1);

t_before = settings.TT.blocktimes.setup;
t_after = settings.TT.blocktimes.afterIC;

% Create empties for all tracks
for tr = 1:Ntracks
    lab = ['track' int2str(tr)];
    timetable.(lab).cancelled(1) = 0;
    timetable.(lab).adjusted_thistrack(1) = 0;
    timetable.(lab).adjusted_start(1) = 0;
    timetable.(lab).adjusted_arrival(1) = 0;
    timetable.(lab).adjusted_departure(1) = 0;
    timetable.(lab).adjusted_finish(1) = 0;
end
    

for tt = 1:Ntrains
    % Which timetable has to be updated?
    if ~traininfo(tt).cancelled
        track = traininfo(tt).newtrack;
        lab = ['track' int2str(track)];
    else
        track = traininfo(tt).track;
        lab = {};
        for tr = 1:Ntracks
            lab(tr) = {['track' int2str(settings.tracks.nr(tr))]};
        end
    end
    
    % Which events come into play?
    ev = traininfo(tt).ev;
    if traininfo(tt).cancelled
        for ll = 1:length(lab)
            timetable.(lab{ll}).cancelled(ev) = 1;
        end
    else
        % Not cancelled, nor deviated. Take the right track!
        timetable.(lab).adjusted_thistrack(ev) = track;
        timetable.(lab).adjusted_arrival(ev(1)) = traininfo(tt).adjusted_entry;
        for ee = 1:length(ev)-1
            timetable.(lab).adjusted_start(ev(ee)) = timetable.(lab).adjusted_arrival(ev(ee))...
                        - timetable.(lab).approachtime(ev(ee)) - t_before;
            timetable.(lab).adjusted_departure(ev(ee)) = timetable.(lab).adjusted_arrival(ev(ee))...
                        + timetable.(lab).running(ev(ee)); 
            timetable.(lab).adjusted_finish(ev(ee)) = timetable.(lab).adjusted_departure(ev(ee))...
                        + t_after;   
            timetable.(lab).adjusted_arrival(ev(ee+1)) = timetable.(lab).adjusted_departure(ev(ee));
        end
        timetable.(lab).adjusted_start(ev(end)) = timetable.(lab).adjusted_arrival(ev(end))...
                    - timetable.(lab).approachtime(ev(end)) - t_before;
        timetable.(lab).adjusted_departure(ev(end)) = timetable.(lab).adjusted_arrival(ev(end))...
                    + timetable.(lab).running(ev(end)); 
        timetable.(lab).adjusted_finish(ev(end)) = timetable.(lab).adjusted_departure(ev(end))...
                    + t_after; 
    end
end

new_timetable = timetable;

end