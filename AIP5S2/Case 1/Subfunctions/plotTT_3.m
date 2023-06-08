function [line, blocks] = plotTT_3(timetable, blocksections, settings, type, include_blocks, firstTime)

% Keep only the rows with '03' or '13' in the 'direction' column
direction_filter = (timetable.direction == 03) | (timetable.direction == 13);
timetable = timetable(direction_filter, :);

% Keep only the rows with 'S1', 'S2', 'C' or 'D' in the 'section' column
section_filter = ismember(blocksections.type, {'S1', 'S2', 'C', 'D'});
blocksections = blocksections(section_filter, :);

try 
    blocksections.closed(1);
catch
    blocksections.closed(1) = 0;
end

if nargin < 6
    firstTime = 0;
end

line = [];
ll = 0;
blocks = [];
bb = 0;

% Generate the structure with lines and blocks, based on the timetable we
% have.

% Run over the complete timetable!
for ee = 1:size(timetable,1)
    
    switch timetable.train_type{ee}
        case 'R'
            C = 'blue';
        case 'IC'
            C = 'black';
        case 'THA'
            C = 'red';
    end
    switch timetable.direction(ee)
        case 03
            L = ':';
        case 12
            L = '--';
        case 13
            L = '-';
    end

    block_id = timetable.blocksection(ee);
    block = find(blocksections.id == block_id);
    % First, draw the line towards the previous value if it is the same
    % train line.
    try 
        if timetable.train_id(ee) == timetable.train_id(ee-1)
            % Connect the arrival at this one and the departure from the
            % previous one. Mainly to detect errors!
            ll = ll+1;
            if timetable.direction(ee) == 13
                line(ll).X1 = blocksections.distance(block) ...
                    - blocksections.length(block);
            else
                line(ll).X1 = blocksections.distance(block);
            end
            line(ll).X2 = line(ll).X1;
            line(ll).Y1 = timetable.departure(ee-1);
            line(ll).Y2 = timetable.arrival(ee);  
            line(ll).Color = C;
            line(ll).LineStyle = L;
        end
    end
    
    % Now the event on this section
    ll = ll+1;
    if timetable.direction(ee) == 13
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
    line(ll).LineStyle = L;
    
    % Also the block
    bb = bb+1;
    if timetable.direction(ee) == 13
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

%% Actual plotting
figure;
hold on

% First, highlight the closed area
switch type
    case 'base'
        if max(timetable.finish) < 1200
            max_time = 1200;
        elseif max(timetable.finish) < 1800
            max_time = 1800;
        elseif max(timetable.finish) < 2400
            max_time = 2400;
        else
            max_time = 3600;
        end
    case 'hour'
        min_time = min([line.Y1, line.Y2]);
        max_time = min_time + 3600;
    case 'complete'
        max_time = settings.disruption.duration*3600;
end

for bb = 1:size(blocksections,1)
    if blocksections.closed(bb)
        x1 = blocksections.distance(bb) - blocksections.length(bb);
        rectangle('Position', [x1 -max_time blocksections.length(bb) max_time*2],...
                    'FaceColor', 0.90*[1 1 1],'EdgeColor','none');
    end
end



if include_blocks
    % Draw rectangles
    % rectangle('Position', [x_left_bottom, y_left_bottom, widht, height])
    for bb = 1:size(blocks,2)
        width = blocks(bb).right - blocks(bb).left;
        height = blocks(bb).top - blocks(bb).bottom;
        p = patch([blocks(bb).left  blocks(bb).right blocks(bb).right blocks(bb).left ],...  
            [blocks(bb).bottom blocks(bb).bottom blocks(bb).top blocks(bb).top],blocks(bb).Color);
        p.FaceAlpha = 0.3;
        p.EdgeAlpha = 0;
        
%         r = rectangle('Position',[blocks(bb).left, blocks(bb).bottom, width, height],...
%             'EdgeColor',blocks(bb).Color,'FaceColor',blocks(bb).Color);
        
    end
end

% Now plot everything
for ll = 1:size(line,2)
    plot([line(ll).X1, line(ll).X2], [line(ll).Y1, line(ll).Y2],'Color',line(ll).Color,'LineStyle',line(ll).LineStyle);
end

% Modify the axes
set(gca,'Ydir','reverse')
set(gca,'xaxislocation','top');
switch type
    case 'base'
        if max(timetable.finish) < 1200
            max_time = 1200;
        elseif max(timetable.finish) < 1800
            max_time = 1800;
        elseif max(timetable.finish) < 2400
            max_time = 2400;
        else
            max_time = 3600;
        end
    case 'hour'
        min_time = min([line.Y1, line.Y2]);
        max_time = min_time + 3600;
    case 'complete'
        max_time = settings.disruption.duration*3600;
end

try
    min_time;
catch
    min_time = 0;
end


xlim([0 max(blocksections.distance)]);
%xlim([0 max(blocksections.distance(strcmp(blocksections.type,'S2')))]);
ylim([min_time max_time]);

% Get interesting axes!
ax = gca;
% set(ax,'XTick',[0, blocksections.distance'],'YGrid','on')
% Where to put the YTicks? => every 15min?
minT = ceil((min_time + firstTime)/900) * 900;
YTickStart = mod(min_time + firstTime,900) + min_time;
maxT = floor(max_time/900) * 900;
maxT = floor((max_time + firstTime)/900) * 900;
YTickEnd = floor((YTickStart + max_time - min_time)/900) * 900; % + firstTime;
YTick = (YTickStart:900:YTickEnd);
YTick_HHMMSS = {};
YTickLab_num = minT:900:maxT;
% First time
startT = firstTime + min_time;
if startT ~= YTickLab_num(1)
    text = timeHHMMSS(startT);
    YTick_HHMMSS{1} = text(1:end-3);
%     YTick =[YTick(1) + startT - YTickLab_num(1), YTick];
    YTick =[min_time, YTick];
    for yy = 2:length(YTick)
        text = timeHHMMSS(YTickLab_num(yy-1));
        YTick_HHMMSS{yy} = text(1:end-3);
    end
else
    for yy = 1:length(YTick)
        text = timeHHMMSS(YTickLab_num(yy));
        YTick_HHMMSS{yy} = text(1:end-3);
    end
end
% Last time
endT = startT + max_time - min_time;
if endT ~= YTick(end)
    text = timeHHMMSS(endT);
    YTick_HHMMSS{end+1} = text(1:end-3);
    YTick =[YTick, max_time + YTick(1) - min_time];
end
% YTick = YTick - 3600;


while YTick(end) == YTick(end-1)
    YTick(end) = [];
    YTick_HHMMSS(end) = [];
end

set(ax,'YTick',YTick,'YTickLabel',YTick_HHMMSS,'YGrid','on')
% set(ax, 'Ydir', 'reverse')
ylabel('Time [H:MM]');
xlabel('Distance [m]');
createTitle = [settings.general.caseName ' ' settings.general.subName ...
                    ' original'];
title(createTitle);


% Show the closed area
try
    blocksections.closed(1);
catch
    blocksections.closed(1) = 0;
end

hold off;



end