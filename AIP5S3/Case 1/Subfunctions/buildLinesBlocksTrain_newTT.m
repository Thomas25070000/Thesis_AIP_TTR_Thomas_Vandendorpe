function [lines, blocks, dots] = buildLinesBlocksTrain_newTT(blocksections, fulltimetable, train_id, track, origtrack, lines, blocks, dots, plotSettings); 

ll = size(lines,2);
bb = size(blocks,2);
dd = size(dots,2);


lab = ['track' int2str(track)];
try
    TTtrack = struct2table(fulltimetable.(lab));
catch
    TTtrack = fulltimetable.(lab);
end
timetable = TTtrack(find([TTtrack.train_id] == train_id),:);


for ee = 1:size(timetable,1)
    
    dontPlot = 0;
    if timetable.cancelled(ee)
        C = 'green';
        style = ':';
    elseif timetable.adjusted_thistrack(ee) ~= track && ~plotSettings.grouptracks
        % Plot the original path if it is moved to another track
        C = 'cyan';
        style = ':';
        if ~plotSettings.origTrack
            dontPlot = 1;
        end
    else
        switch timetable.direction(ee)
            case {02, 12}
                C = 'red';
            case {03, 13}
                C = 'blue';
        end
        style = '-';
    end
    if timetable.adjusted_thistrack(ee)== 1
        C = 'black';
        style = '-';
    end
   
    if ~dontPlot
        block = timetable.blocksection(ee);

        % First, draw the line towards the previous value if it is the same
        % train line.
        try 
            if timetable.train_id(ee) == timetable.train_id(ee-1)
                % Connect the arrival at this one and the departure from the
                % previous one. Mainly to detect errors!
                ll = ll+1;
                if timetable.direction(ee) == 12 || timetable.direction(ee) == 13
                    lines(ll).X1 = blocksections.distance(block) ...
                        - blocksections.length(block);
                else
                    lines(ll).X1 = blocksections.distance(block);
                end
                lines(ll).X2 = lines(ll).X1;
                if timetable.cancelled(ee)...
                        || (timetable.adjusted_thistrack(ee) ~= track && ~plotSettings.grouptracks)
                    lines(ll).Y1 = timetable.departure(ee-1);
                    lines(ll).Y2 = timetable.arrival(ee);  
                else
                    lines(ll).Y1 = timetable.adjusted_departure(ee-1);
                    lines(ll).Y2 = timetable.adjusted_arrival(ee);  
                end
                lines(ll).Color = C;
                lines(ll).Style = style;
            elseif ~(timetable.adjusted_thistrack(ee) ~= track && ~plotSettings.grouptracks)
                dd = dd+1;
                if timetable.direction(ee) == 12 || timetable.direction(ee) == 13
                    dots(dd).X = blocksections.distance(block) ...
                        - blocksections.length(block);
                else
                    dots(dd).X = blocksections.distance(block);
                end
                dots(dd).Y = timetable.arrival(ee);
                dots(dd).Color = C;
                if timetable.adjusted_thistrack(ee) == track && timetable.thistrack(ee)
                    dots(dd).Fill = 1;
                else
                    dots(dd).Fill = 1;
                end
            end
        catch
            dd = dd+1;
            if timetable.direction(ee) == 12 || timetable.direction(ee) == 13
                dots(dd).X = blocksections.distance(block) ...
                    - blocksections.length(block);
            else
                dots(dd).X = blocksections.distance(block);
            end
            dots(dd).Y = timetable.arrival(ee);
            dots(dd).Color = C;
            if timetable.adjusted_thistrack(ee) == track && timetable.thistrack(ee)
                dots(dd).Fill = 1;
            else
                dots(dd).Fill = 1;
            end
        end

        % Now the event on this section
        ll = ll+1;
        if timetable.direction(ee) == 12 || timetable.direction(ee) == 13
            lines(ll).X1 = blocksections.distance(block) ...
                        - blocksections.length(block);
            lines(ll).X2 = blocksections.distance(block);
        else
            lines(ll).X1 = blocksections.distance(block);
            lines(ll).X2 = blocksections.distance(block) ...
                        - blocksections.length(block);
        end
        if timetable.cancelled(ee)...
                    || (timetable.adjusted_thistrack(ee) ~= track && ~plotSettings.grouptracks)
            lines(ll).Y1 = timetable.arrival(ee);
            lines(ll).Y2 = timetable.departure(ee);
        else
            lines(ll).Y1 = timetable.adjusted_arrival(ee);
            lines(ll).Y2 = timetable.adjusted_departure(ee);
        end
        lines(ll).Color = C;
        lines(ll).Style = style;

        flag_add_block = 1;
        if timetable.cancelled(ee)...
                || (timetable.adjusted_thistrack(ee) ~= track && ~plotSettings.grouptracks)
            flag_add_block = 0;
        end

        if flag_add_block
            bb = bb+1;
            if timetable.direction(ee) == 12 || timetable.direction(ee) == 13
                blocks(bb).left = lines(ll).X1;
                blocks(bb).right = lines(ll).X2;
            else
                blocks(bb).left = lines(ll).X2;
                blocks(bb).right = lines(ll).X1;
            end
            blocks(bb).bottom = timetable.adjusted_start(ee);
            blocks(bb).top = timetable.adjusted_finish(ee);
            blocks(bb).Color = C;
        end
    end
    
end


end