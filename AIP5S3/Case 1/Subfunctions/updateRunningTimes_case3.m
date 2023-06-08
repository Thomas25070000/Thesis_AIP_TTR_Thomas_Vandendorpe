function timetable = updateRunningTimes_case3(timetable, runningtimes, settings, traininfo)

Ntracks = size(settings.tracks,1);
for tr = 1:Ntracks
    if settings.tracks.vreg(tr) ~= settings.tracks.vdis(tr) && ~settings.tracks.closed(tr)
        % We have to update the timetable!
        labtrack = ['track' int2str(settings.tracks.nr(tr))];
        for tt = 1:size(traininfo,2)
            events = traininfo(tt).ev;
            dir = traininfo(tt).dir;
            type = traininfo(tt).type;
            switch type
                case {'IC', 'THA'}
                    type = 'IC';
                case {'R','L','P'}
                    type = 'L';
                otherwise
                    error('Unrecognized train type');
            end
            labtrain = [type int2str(dir)];
            timetable.(labtrack).running(events) = runningtimes.(labtrain).disrupted(tr,:);
        end
    end
end
    
    
end