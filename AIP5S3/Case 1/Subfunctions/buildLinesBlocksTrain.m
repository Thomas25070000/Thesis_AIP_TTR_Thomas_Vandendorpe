function [lines_train, blocks_train] = buildLinesBlocksTrain(blocksections, fulltimetable, train_id, track)

line = [];
ll = 0;
blocks = [];
bb = 0;

lab = ['track' int2str(track)];
try
    TTtrack = struct2table(fulltimetable.(lab));
catch
    TTtrack = fulltimetable.(lab);
end
timetable = TTtrack(find([TTtrack.train_id] == train_id),:);


for ee = 1:size(timetable,1)
    
    switch timetable.train_type{ee}
        case 'R'
            C = 'blue';
        case 'IC'
            C = 'black';
        case 'THA'
            C = 'red';
    end
    block = timetable.blocksection(ee);
    
    % First, draw the line towards the previous value if it is the same
    % train line.
    try 
        if timetable.train_id(ee) == timetable.train_id(ee-1)
            % Connect the arrival at this one and the departure from the
            % previous one. Mainly to detect errors!
            ll = ll+1;
            if timetable.direction(ee) == 1
                line(ll).X1 = blocksections.distance(block) ...
                    - blocksections.length(block);
            else
                line(ll).X1 = blocksections.distance(block);
            end
            line(ll).X2 = line(ll).X1;
            line(ll).Y1 = timetable.departure(ee-1);
            line(ll).Y2 = timetable.arrival(ee);  
            line(ll).Color = C;
        end
    end
    
    % Now the event on this section
    ll = ll+1;
    if timetable.direction(ee) == 1
        line(ll).X1 = blocksections.distance(block) ...
                    - blocksections.length(block);
        line(ll).X2 = blocksections.distance(block);
    else
        line(ll).X1 = blocksections.distance(block);
        line(ll).X2 = blocksections.distance(block) ...
                    - blocksections.length(block);
    end
    line(ll).Y1 = timetable.arrival(ee);
    line(ll).Y2 = timetable.departure(ee);
    line(ll).Color = C;
    
    % Also the block
    bb = bb+1;
    if timetable.direction(ee) == 1
        blocks(bb).left = line(ll).X1;
        blocks(bb).right = line(ll).X2;
    else
        blocks(bb).left = line(ll).X2;
        blocks(bb).right = line(ll).X1;
    end
    blocks(bb).bottom = timetable.start(ee);
    blocks(bb).top = timetable.finish(ee);
    blocks(bb).Color = C;
end

lines_train = line;
blocks_train = blocks;



end