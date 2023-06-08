function timetable = updateRunningTimesAndTT_case3(timetable, runningtimes, settings, traininfo, blocksections)

Ntracks = size(settings.tracks,1);
t_setup = settings.TT.blocktimes.setup;
if settings.TT.blocktimes.afterR ~= settings.TT.blocktimes.afterIC
    error('Issue with release times: not equal for both train types');
end
t_after = settings.TT.blocktimes.afterR;


for tr = 1:Ntracks
    labtrack = ['track' int2str(settings.tracks.nr(tr))];
    timetable.(labtrack).approachtime(1) = 0;
    
%     if settings.tracks.vreg(tr) ~= settings.tracks.vdis(tr) && ~settings.tracks.closed(tr)
    if ~settings.tracks.closed(tr)
        % We have to update the timetable!
        
        for tt = 1:size(traininfo,2)
            events = traininfo(tt).ev;
            dir = traininfo(tt).dir;
            type = traininfo(tt).type;
            switch type
                case {'IC', 'THA'}
                    type = 'IC';
                case {'R','L','P','S'}
                    type = 'L';
                otherwise
                    error('Unrecognized train type');
            end
            labtrain = [type int2str(dir)];
            if dir == 2
                labtrain = 'IC02';
            end
            if dir == 3
                labtrain = 'IC03';
            end
            if dir == 2
                indices = [5 6 7 8 9 10 11];
            end
            if dir == 3
                indices = [1 2 3 8 9 10 11];
            end
            if dir == 12
                indices = [1 2 3 4 5 6 7];
            end
            if dir == 13
                indices = [11 10 9 8 3 2 1];
            end
             
            timetable.(labtrack).running(events) = runningtimes.(labtrain).disrupted(tr,indices);
            
            % Fill the approach times and update the timetable, including
            % blocking times.
            [app, ~] = returnFirstApproachTime_case3(blocksections, settings, type, dir, tr);
            
            timetable.(labtrack).start(events(1)) = timetable.(labtrack).arrival(events(1)) - app - t_setup;
            timetable.(labtrack).approachtime(events(1)) = app;
            
            for ee = 1:length(events)-1
                timetable.(labtrack).departure(events(ee)) = timetable.(labtrack).arrival(events(ee))...
                                + timetable.(labtrack).running(events(ee));
                timetable.(labtrack).finish(events(ee)) = timetable.(labtrack).departure(events(ee)) + t_after;
                timetable.(labtrack).arrival(events(ee+1)) = timetable.(labtrack).departure(events(ee));        
                app = timetable.(labtrack).running(events(ee));
                timetable.(labtrack).approachtime(events(ee+1)) = app;
                timetable.(labtrack).start(events(ee+1)) = timetable.(labtrack).arrival(events(ee+1)) ...
                                - app - t_setup;
            end
            timetable.(labtrack).departure(events(end)) = timetable.(labtrack).arrival(events(end))...
                                + timetable.(labtrack).running(events(end));
            timetable.(labtrack).finish(events(end)) = timetable.(labtrack).departure(events(end)) + t_after;                

        end

    end
    
end




    
    
end