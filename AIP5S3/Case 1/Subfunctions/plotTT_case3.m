function [line, blocks] = plotTT_case3(timetable, blocksections, traininfo, settings, plotSettings)

try 
    blocksections.closed(1);
catch
    blocksections.closed(1) = 0;
end

try
    firstTime = plotSettings.firstTime;
catch
    firstTime = 0;
end
try
    include_blocks = plotSettings.include_blocks;
catch
    include_blocks = 1;
end

lines = [];
ll = 0;
blocks = [];
bb = 0;

% Remove abundant train info
trainsTT = unique([timetable.track1.train_id]);
remove = find(~ismember([traininfo.id],trainsTT));
traininfo(remove) = [];

% Generate the structure with lines and blocks, based on the timetable we
% have.
 
% Do we want all on one, or separated?
if plotSettings.grouptracks
    % All on one!
    for tt = 1:size(traininfo,2)
        train_id = traininfo(tt).id;
        track = traininfo(tt).track;
        if ismember(track, plotSettings.tracks)
            [lines_train, blocks_train] = buildLinesBlocksTrain(blocksections, timetable, train_id, track);
            if ~isempty(lines)
                f = fieldnames(lines_train);
                nrLines = size(lines,2);
                for ll = 1:size(lines_train,2)
                    for i = 1:length(f)
                        lines(nrLines+ll).(f{i}) = lines_train(ll).(f{i});
                    end
                end
                
                
            else
                lines = lines_train;
            end
            if ~isempty(blocks) && include_blocks
                f = fieldnames(blocks_train);
                nrBlocks = size(blocks,2);
                for bb = 1:size(blocks_train,2)
                    for i = 1:length(f)
                        blocks(nrBlocks+bb).(f{i}) = blocks_train(bb).(f{i});
                    end
                end
            elseif include_blocks
                blocks = blocks_train;
            end
        end
    end
else
    % Repeat the procedure above for each and every track.
    for tr = 1:length(plotSettings.tracks)
        blocks = [];
        lines = [];

        for tt = 1:size(traininfo,2)
            train_id = traininfo(tt).id;
            track = traininfo(tt).track;
            if ismember(track, plotSettings.tracks(tr))
                [lines_train, blocks_train] = buildLinesBlocksTrain(blocksections, timetable, train_id, track);
                if ~isempty(lines)
                    f = fieldnames(lines_train);
                    nrLines = size(lines,2);
                    for ll = 1:size(lines_train,2)
                        for i = 1:length(f)
                            lines(nrLines+ll).(f{i}) = lines_train(ll).(f{i});
                        end
                    end
                else
                    lines = lines_train;
                end
                if ~isempty(blocks) && include_blocks
                    f = fieldnames(blocks_train);
                    nrBlocks = size(blocks,2);
                    for bb = 1:size(blocks_train,2)
                        for i = 1:length(f)
                            blocks(nrBlocks+bb).(f{i}) = blocks_train(bb).(f{i});
                        end
                    end
                elseif include_blocks
                    blocks = blocks_train;
                end
            end
        end
        
        lab = ['track' int2str(plotSettings.tracks(tr))];
        allblocks.(lab) = blocks;
        alllines.(lab) = lines;
    end
end


%% Actual plotting
if plotSettings.grouptracks
    line = lines;
    blocks = blocks;
    plottitle = [settings.general.subName ' original (tracks '];
    for tr = 1:length(plotSettings.tracks)-1
        plottitle = [plottitle int2str(plotSettings.tracks(tr)) ', '];
    end
    plottitle = [plottitle int2str(plotSettings.tracks(end)) ')'];    
    
    figure;
    hold on

    % First, highlight the closed area
    switch plotSettings.type
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
        plot([line(ll).X1, line(ll).X2], [line(ll).Y1, line(ll).Y2],'Color',line(ll).Color);
    end

    % Modify the axes
    set(gca,'Ydir','reverse')
    set(gca,'xaxislocation','top');
    switch plotSettings.type
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
    while  YTick(end) == YTick(end-1)
        YTick(end) = [];
    end

    set(ax,'YTick',YTick,'YTickLabel',YTick_HHMMSS,'YGrid','on')
    % set(ax, 'Ydir', 'reverse')
    ylabel('Time [H:MM]');
    xlabel('Distance [m]');
    createTitle = [settings.general.caseName ' ' settings.general.subName ...
                        ' original'];
%     title(createTitle);
    title(plottitle);


    % Show the closed area
    try
        blocksections.closed(1);
    catch
        blocksections.closed(1) = 0;
    end

    hold off;

    
else
    % Generate a plot for each and every track!
    for tr = 1:length(plotSettings.tracks)
        track = plotSettings.tracks(tr);
        lab = ['track' int2str(track)];
        
        line = alllines.(lab);
        blocks = allblocks.(lab);
        plottitle = [settings.general.subName ' original (track ' int2str(track) ')'];
        
        figure;
        hold on

        % First, highlight the closed area
        switch plotSettings.type
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
            plot([line(ll).X1, line(ll).X2], [line(ll).Y1, line(ll).Y2],'Color',line(ll).Color);
        end

        % Modify the axes
        set(gca,'Ydir','reverse')
        set(gca,'xaxislocation','top');
        switch plotSettings.type
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
        end

        set(ax,'YTick',YTick,'YTickLabel',YTick_HHMMSS,'YGrid','on')
        % set(ax, 'Ydir', 'reverse')
        ylabel('Time [H:MM]');
        xlabel('Distance [m]');
        createTitle = [settings.general.caseName ' ' settings.general.subName ...
                            ' original'];
    %     title(createTitle);
        title(plottitle);


        % Show the closed area
        try
            blocksections.closed(1);
        catch
            blocksections.closed(1) = 0;
        end

        hold off;
    end
end

end