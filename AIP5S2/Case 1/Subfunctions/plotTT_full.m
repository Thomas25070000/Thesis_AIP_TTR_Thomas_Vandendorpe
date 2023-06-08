function [line_new, blocks_new] = plotTT_full(new_timetable, blocksections, settings, type, include_blocks, directions, firstTime)

% % Keep only the rows with '02' or '12' in the 'direction' column
% direction_filter = (new_timetable.direction == 02) | (new_timetable.direction == 12);
% new_timetable = new_timetable(direction_filter, :);
% 
% % Keep only the rows with 'S1', 'S2', 'B' or 'D' in the 'section' column
% section_filter = ismember(blocksections.type, {'S1', 'S2', 'B', 'D'});
% blocksections = blocksections(section_filter, :);

%new_timetable(new_timetable.direction == 13, :) = [];
try
    noDots = settings.noDots;
catch
    noDots = 0;
end

if nargin < 7
    try
		firstTime = settings.firstTime;
	catch
		firstTime = 0;
	end	
end

timetable = new_timetable;

try 
    blocksections.closed(1);
catch
    blocksections.closed(1) = 0;
end

try
    directions;
catch
    directions = [02 03 12 13];
end


line = [];
ll = 0;
blocks = [];
bb = 0;
dots = [];      % For the original, planned arrivals.
dd = 0; 

% Generate the structure with lines and blocks, based on the timetable we
% have.

% Run over the complete timetable!
for ee = 1:size(timetable,1)
    

    switch timetable.direction(ee)
        case {2,12}
            C = 'r';
        case {3}
            C = 'b';
        case {13}
            C = 'black';
        case 'THA'
            C = 'red';
    end
    
    switch timetable.direction(ee)
        case {2, 12}
            color = 'r';
        case {3}
            color = 'b';
        case {13}
            color = 'black';
    end

%         % Only show specific direction completely
%         if ismember(timetable.direction(ee), directions)
%             style = '-';
%         else
%             style = ':';
%         end
    if timetable.cancelled(ee)
        C = 'green';
        color = 'green';
        style = ':';
    elseif settings.deviation.available && timetable.deviated(ee)
        C = 'green';
        style = ':';
    else
    style = '-';
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
            if timetable.direction(ee) == 12 | timetable.direction(ee) == 13
                line(ll).X1 = blocksections.distance(block) ...
                    - blocksections.length(block);
            else
                line(ll).X1 = blocksections.distance(block);
            end
            line(ll).X2 = line(ll).X1;
            if timetable.cancelled(ee) || (settings.deviation.available && timetable.deviated(ee))
                line(ll).Y1 = timetable.departure(ee-1);
                line(ll).Y2 = timetable.arrival(ee);  
            else
                line(ll).Y1 = timetable.adjusted_departure(ee-1);
                line(ll).Y2 = timetable.adjusted_arrival(ee);  
            end
            line(ll).Color = C;
            line(ll).Style = style;
            line(ll).Color = color;
        else
            dd = dd+1;
            if timetable.direction(ee) == 12 | timetable.direction(ee) == 13
                dots(dd).X = blocksections.distance(block) ...
                    - blocksections.length(block);
                dots(dd).Y = timetable.arrival(ee);
            else
                %dots(dd).X = blocksections.distance(block);
                dots(dd).X = 4400;
                dots(dd).Y = timetable.arrival(ee+4);
            end
            
            dots(dd).Color = C;
        end
    catch
        dd = dd+1;
        if timetable.direction(ee) == 12 | timetable.direction(ee) == 13
            dots(dd).X = blocksections.distance(block) ...
                - blocksections.length(block);
            dots(dd).Y = timetable.arrival(ee);
        else
            %dots(dd).X = blocksections.distance(block);
            dots(dd).X = 4400;
            dots(dd).Y = timetable.arrival(ee+4);
        end
        
        dots(dd).Color = C;
    end
    
    % Now the event on this section
    ll = ll+1;
    if timetable.direction(ee) == 12 | timetable.direction(ee) == 13
        line(ll).X1 = blocksections.distance(block) ...
                    - blocksections.length(block);
        line(ll).X2 = blocksections.distance(block);
    else
        line(ll).X1 = blocksections.distance(block);
        line(ll).X2 = blocksections.distance(block) ...
                    - blocksections.length(block);
    end
    
    if timetable.cancelled(ee) || (settings.deviation.available && timetable.deviated(ee))
        line(ll).Y1 = timetable.arrival(ee);
        line(ll).Y2 = timetable.departure(ee);
    else
        line(ll).Y1 = timetable.adjusted_arrival(ee);
        line(ll).Y2 = timetable.adjusted_departure(ee);
    end
    line(ll).Color = C;
    line(ll).Style = style;
    line(ll).Color = color;
    
    % Also the block
    include_block = 0;
    if ismember(timetable.direction(ee),directions)
        include_block = 1;
    elseif blocksections.closed(timetable.blocksection(ee))
        include_block = 1;
    end
        
    
    % Flag to satisty conditions of adding a block or not
    flag_add_block = include_block;
    if timetable.cancelled(ee) || (settings.deviation.available && timetable.deviated(ee))
        flag_add_block = 0;
%     elseif floor(timetable.train_id(ee)/1000) > settings.disruption.duration
    elseif timetable.adjusted_start(ee) > settings.disruption.duration*3600
        if ~ismember(timetable.direction(ee),directions)
            flag_add_block = 0;
        end
    end
        
        
    if flag_add_block
        bb = bb+1;
        if timetable.direction(ee) == 12 | timetable.direction(ee) == 13
            blocks(bb).left = line(ll).X1;
            blocks(bb).right = line(ll).X2;
        else
            blocks(bb).left = line(ll).X2;
            blocks(bb).right = line(ll).X1;
        end
        blocks(bb).bottom = timetable.adjusted_start(ee);
        blocks(bb).top = timetable.adjusted_finish(ee);
        blocks(bb).Color = C;
    end
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
        max_time = (settings.disruption.duration + settings.disruption.nr_hours_after) * 3600;
end

for bb = 1:size(blocksections,1)
    if blocksections.closed(bb)
        x1 = blocksections.distance(bb) - blocksections.length(bb);
        rectangle('Position', [x1 0 blocksections.length(bb) ...
            (settings.disruption.duration+settings.disruption.nr_hours_after)*3600],'FaceColor', 0.90*[1 1 1],'EdgeColor','none');
    end
end


if include_blocks
    % Draw rectangles
    % rectangle('Position', [x_left_bottom, y_left_bottom, widht, height])
    for bb = 1:size(blocks,2)
        width = blocks(bb).right - blocks(bb).left;
        height = blocks(bb).top - blocks(bb).bottom;
        
        if width <=0 | height <=0
            disp('Problem');
        end
		p = patch([blocks(bb).left  blocks(bb).right blocks(bb).right blocks(bb).left ],...  
            [blocks(bb).bottom blocks(bb).bottom blocks(bb).top blocks(bb).top],blocks(bb).Color);
        p.FaceAlpha = 0.3;
        p.EdgeAlpha = 0;
        
  %      rectangle('Position',[blocks(bb).left, blocks(bb).bottom, width, height],'EdgeColor',blocks(bb).Color);
    end
end

% Now plot everything
for ll = 1:size(line,2)
    plot([line(ll).X1, line(ll).X2], [line(ll).Y1, line(ll).Y2],'Color',line(ll).Color, 'LineStyle',line(ll).Style);
end

if ~noDots
    for dd = 1:size(dots,2)
        scatter(dots(dd).X,dots(dd).Y,20,dots(dd).Color,'filled');
    end
end

% Modify the axes
set(gca,'Ydir','reverse')
set(gca,'xaxislocation','top');


try
    min_time;
catch
    min_time = 0;
end

xlim([0 max(blocksections.distance)]);
%xlim([0 max(blocksections.distance(strcmp(blocksections.type,'S2')))]);
%ylim([min_time max_time]);
%ylim([min_time min_time+21600]);
ylim([0 7200]);

% Plot horizontal lines at the hour borders
% hours = floor(max_time/3600);
% if ~noDots
%     for hh = 1:hours
%         plot([0 max(blocksections.distance)],[hh*3600 hh*3600],':r');
%     end
% end

% Get interesting axes!
ax = gca;


% Set YTick values
YTickStart = 0;
YTickEnd = 21600;
YTickStep = 1800;
YTick = YTickStart:YTickStep:YTickEnd;

% Set YTickLabel values in HH:MM format
YTick_HHMMSS = {};
for t = YTick/YTickStep
    hours = floor((t+1)/2); % Each interval of 1800 seconds (30 min) represents 1 hour
    minutes = mod(t+1, 2) * 30; % The remainder of the division by 2 indicates if it is the first or second half of the hour (0 or 30 minutes)
    YTick_HHMMSS{end+1} = sprintf('%02d:%02d', hours+15, minutes); % Add 15 to hours to start at 15:30
end

set(ax,'YTick',YTick,'YTickLabel',YTick_HHMMSS,'YGrid','on')
% set(ax, 'Ydir', 'reverse')
ylabel('Time [H:MM]');
xlabel('Distance [m]');
% createTitle = [settings.general.caseName ' ' settings.general.subName ...
%                     ' adjusted (v' int2str(settings.general.modelVersion) ')'];
createTitle = [settings.general.caseName ' ' settings.general.subName ...
                    ' adjusted'];
title('Adjusted timetable during first 2 hours using a FIFO heuristic, wCancel=3600.');



% Show the closed area
try
    blocksections.closed(1);
catch
    blocksections.closed(1) = 0;
end

hold off;



line_new = line;
blocks_new = block;

end