function [measures, statistics] = deriveMeasuresFIFO(timetable, cancelled, settings, orig_delays, extra_delays)

saveStats = settings.saveStats;


% Already give some of the statistics here!
% Generate a delay table
delays = [];
row = 0;
trains = unique(timetable.train_id);
for tt = 1:length(trains)
    train_ev = timetable(find(timetable.train_id == trains(tt)),:);
    deviated(tt) = train_ev.deviated(1);
    temp = timetable.arrival(find(timetable.train_id == trains(tt)));
    t_closed(tt) = temp(1);
end

direction = mod(trains,10);
for dd = 1:length(orig_delays)
    % Also generate HH:MM:SS
    if ~cancelled(dd)
        time = orig_delays(dd);
        extra_time = extra_delays(dd);
        row = row + 1;
        delays(row).train_id = trains(dd);
        delays(row).orig_delay = time;
        delays(row).exta_delay = time;
        delays(row).orig_delayHHMMSS = timeHHMMSS(time);
    end
end

w_dir0 = settings.weights.pax_dir0 / settings.weights.pax_dir1;
w_dir1 = 1;

if w_dir0 ~= w_dir1
    Objective = 0;
    trains = unique(timetable.train_id);
    
    for tt = 1:length(trains)
        if cancelled(tt)
            if mod(trains(tt),10) == 1
                Objective = Objective + settings.weights.cancel * w_dir1;
            else
                Objective = Objective + settings.weights.cancel * w_dir0;
            end
        else
            if mod(trains(tt),10) == 1
                Objective = Objective + extra_delays(tt) * w_dir1;
            else
                Objective = Objective + extra_delays(tt) * w_dir0;
            end
        end
    end
end





delays_nocancel = [];

if ~isempty(delays)
    delays = struct2table(delays);
    statistics.delays = delays;

    delays_nocancel = orig_delays(find(~cancelled & ~deviated'));
    extra_delays_nocancel = extra_delays(find(~cancelled & ~deviated'));
    trains_nocancel = trains(find(~cancelled & ~deviated'));
    direction_nocancel = direction(find(~cancelled & ~deviated'));
end

if ~isempty(delays_nocancel)

    max_delay_rows = find(delays_nocancel == max(delays_nocancel));
    max_delay_row = max_delay_rows(1);
    delays_dir0 = delays_nocancel(find(direction_nocancel == 0));
    max_delay_rows_dir0 = find(delays_dir0 == max(delays_dir0));
    max_delay_row_dir0 = max_delay_rows_dir0(1);
    delays_dir1 = delays_nocancel(find(direction_nocancel == 1));
    max_delay_rows_dir1 = find(delays_dir1 == max(delays_dir1));
    max_delay_row_dir1 = max_delay_rows_dir1(1);

    % Max. delay?
    statistics.maxDelay = delays_nocancel(max_delay_row);
    statistics.maxDelay_HHMMSS = timeHHMMSS(delays_nocancel(max_delay_row));
    statistics.maxDelay_train = trains_nocancel(max_delay_rows);
    statistics.maxDelay_0 = delays_dir0(max_delay_row_dir0);
    statistics.maxDelay_HHMMSS_0 = timeHHMMSS(delays_dir0(max_delay_row_dir0));
    statistics.maxDelay_1 = delays_dir1(max_delay_row_dir1);
    statistics.maxDelay_HHMMSS_1 = timeHHMMSS(delays_dir1(max_delay_row_dir1));

    % Total delay
    total_delay = sum(delays_nocancel);
    total_extra_delay = sum(extra_delays_nocancel);
    statistics.totalDelay = total_delay;
    statistics.totalExtraDelay = total_extra_delay;
    statistics.totalDelay_HHMMSS = timeHHMMSS(total_delay);
    statistics.totalExtraDelay_HHMMSS = timeHHMMSS(total_extra_delay);
    total_delay_0 = sum(delays_dir0);
    statistics.totalDelay_0 = sum(delays_dir0);
    statistics.totalDelay_HHMMSS_0 = timeHHMMSS(sum(delays_dir0));
    total_delay_1 = sum(delays_dir1);
    statistics.totalDelay_1 = sum(delays_dir1);
    statistics.totalDelay_HHMMSS_1 = timeHHMMSS(sum(delays_dir1));

    % Average delay
    nr_delay = sum(delays_nocancel > 0);
    av_delay = total_delay / nr_delay;
    statistics.averageDelay = av_delay;
    statistics.averageDelay_HHMMSS = timeHHMMSS(av_delay);
    nr_delay_0 = sum(delays_dir0 > 0);
    av_delay_0 = total_delay_0 / nr_delay_0;
    statistics.averageDelay_0 = av_delay_0;
    statistics.averageDelay_HHMMSS_0 = timeHHMMSS(av_delay_0);
    nr_delay_1 = sum(delays_dir1 > 0);
    av_delay_1 = total_delay_1 / nr_delay_1;
    statistics.averageDelay_1 = av_delay_1;
    statistics.averageDelay_HHMMSS_1 = timeHHMMSS(av_delay_1);
else
    statistics.maxDelay = 0;
    statistics.maxDelay_HHMMSS = 0;
    statistics.maxDelay_0 = 0;
    statistics.maxDelay_HHMMSS_0 = 0;
    statistics.maxDelay_1 = 0;
    statistics.maxDelay_HHMMSS_1 = 0;
    statistics.maxDelay_train = 0;
    statistics.totalDelay = 0;
    statistics.totalDelay_HHMMSS = 0;
    statistics.totalDelay_0 = 0;
    statistics.totalDelay_HHMMSS_0 = 0;
    statistics.totalDelay_1 = 0;
    statistics.totalDelay_HHMMSS_1 = 0;
    statistics.averageDelay = 0;
    statistics.averageDelay_HHMMSS = 0;
    statistics.averageDelay_0 = 0;
    statistics.averageDelay_HHMMSS_0 = 0;
    statistics.averageDelay_1 = 0;
    statistics.averageDelay_HHMMSS_1 = 0;
end


folder = [pwd '\Output'];

if saveStats
    try
%         filename = [settings.general.caseName ' ' settings.general.subName ...
%                 ' (v' int2str(settings.general.modelVersion) ') - measures.txt'];
        filename = [settings.general.caseName ' ' settings.general.subName ...
                  ' - measures.txt'];
    catch
        filename = generateFilename();
    end
    folderfile = [folder '\' filename];

    fid=fopen(folderfile,'wt');
    
    if fid <= 0
        % Problem with opening the file! Mostly because there is no folder
        % present.
        mkdir('Output');
        fid=fopen(folderfile,'wt');
    end
end

trains = unique(timetable.train_id);

measures.cancelledtrains = trains(find(cancelled));


statistics.cancelledTrains = trains(find(cancelled));
statistics.nrCancelled = sum(cancelled);

measures.deviatedtrains = trains(find(deviated));
statistics.deviatedTrains = trains(find(deviated));
statistics.nrDeviated = sum(deviated);



fprintf('***Measures***\n\n');

if saveStats
    fprintf(fid,'***Measures***\n\n');
end

%% First: cancellations
cancelled_dir0 = cancelled(direction == 0);
cancelled_dir1 = cancelled(direction == 1);
text0 = [int2str(sum(cancelled)) ' trains have to be cancelled\n \t - direction 0: '...
              int2str(sum(cancelled_dir0))];
text1 = ['\n\t - direction 1: ' int2str(sum(cancelled_dir1))];

statistics.nrCancelled_dir0 = sum(cancelled_dir0);
statistics.nrCancelled_dir1 = sum(cancelled_dir1);

text_dir0 = [];
text_dir1 = [];
for cc = 1:length(cancelled)
    if cancelled(cc) == 1
        rows = timetable(find(timetable.train_id == trains(cc)),:);
        type = rows.train_type{1};
        time = rows.arrival(1);
        HHMMSS = [int2str(floor(time/3600)) ':' int2str(floor(mod(time,3600))/60) ':' int2str(mod(time,60))];
        text = ['\n\t\t Train ' type int2str(trains(cc)) ' of ' HHMMSS];
        
        if direction(cc) == 1
            text_dir1 = [text_dir1 text];
        else
            text_dir0 = [text_dir0 text];
        end
    end
end
fprintf([text0 text_dir0 text1 text_dir1]);
if saveStats
    fprintf(fid, [text0 text_dir0 text1 text_dir1]);
end

%% If done: deviations
if settings.deviation.available
    deviated_dir0 = deviated(direction == 0);
    deviated_dir1 = deviated(direction == 1);
    text0 = ['\n\n' int2str(sum(deviated)) ' trains have to be deviated\n \t - direction 0: '...
                  int2str(sum(deviated_dir0))];
    text1 = ['\n\t - direction 1: ' int2str(sum(deviated_dir1))];

    statistics.nrDeviated_dir0 = sum(deviated_dir0);
    statistics.nrDeviated_dir1 = sum(deviated_dir1);

    text_dir0 = [];
    text_dir1 = [];
    
    for cc = 1:length(deviated)
        if deviated(cc) == 1
            rows = timetable(find(timetable.train_id == trains(cc)),:);
            type = rows.train_type{1};
            time = rows.arrival(1);
            HHMMSS = [int2str(floor(time/3600)) ':' int2str(floor(mod(time,3600))/60) ':' int2str(mod(time,60))];
            text = ['\n\t\t Train ' type int2str(trains(cc)) ' of ' HHMMSS];

            if direction(cc) == 1
                text_dir1 = [text_dir1 text];
            else
                text_dir0 = [text_dir0 text];
            end
        end
    end
    
    fprintf([text0 text_dir0 text1 text_dir1]);
    if saveStats
        fprintf(fid, [text0 text_dir0 text1 text_dir1]);
    end
    
    
end

%% Order switching?
text = [];
if settings.constraints.fixEntranceOrder
    text = '\nNo order switching allowed\n';
elseif (isempty(freeOrders_1) && isempty(freeOrders_0))
    text = '\nOrder switching allowed, but not done\n';
else
    text = '\nChanged orders:\n';
    if ~isempty(freeOrders_0)
        text = [text '\t- direction 0: '];
        for oo = 1:size(freeOrders_0,1)
            text = [text '(' int2str(freeOrders_0(oo,1)) ', ' ...
                        int2str(freeOrders_0(oo,2)) '); '];
        end
        text = [text '\n'];
    end
    if ~isempty(freeOrders_1)
        text = [text '\t- direction 1: '];
        for oo = 1:size(freeOrders_1,1)
            text = [text '(' int2str(freeOrders_1(oo,1)) ', ' ...
                        int2str(freeOrders_1(oo,2)) '); '];
        end
        text = [text '\n'];
    end
end
fprintf(text);
if saveStats
    fprintf(fid, text);
end

%% Next: for the remaining trains, what order?
remove = find(cancelled == 1);
trains(remove) = [];
t_closed(remove) = [];
direction(remove) = [];

% Sort the remaining ones according to their time of entry on the closed
% segment.
[t_closed, I] = sortrows(t_closed);
trains_ordered = trains(I);
direction_ordered = direction(I);

if ~isempty(direction_ordered)
    text = ['\n\n Resulting pattern of trains: '];
    text_1 = [];    % Train types
    text_2 = ['\n\t\t\t'];    % Orders
    for tt = 1:length(I)
        text_1 = [text_1 int2str(direction_ordered(tt)) ' '];
        rows = timetable(find(timetable.train_id == trains_ordered(tt)),:);
        type = rows.train_type{1};
        text_2 = [text_2 ' ' type];
    end
    fprintf([text text_1 text_2]);
    if saveStats
        fprintf(fid, [text text_1 text_2]);
    end

    % Frequency pattern?
    text = ['\n\n First direction: ' int2str(direction_ordered(1))];
    fprintf(text);
    if saveStats
        fprintf(fid, text);
    end

    dir = direction_ordered;
    freq = [];
    dd = 2;
    ff = 1;
    hour = 0;
    new_hour = [];
    % while dd <= length(dir)
    %     if dir(dd) == dir(dd-1)
    %         % We are still in the same direction!
    %         % Remove direction.
    %         dir(dd) = [];
    %         ff = ff + 1;
    %         if floor(t_closed(dd)/3600) > hour
    %             % We go to another hour!
    %             hour = hour + 1;
    %             new_hour(length(freq)) = 1;
    %         end
    %         t_closed(dd) = [];            
    %     else
    %         % Change of types!
    %         % Advance to next spot.
    %         dd = dd + 1;
    %         % Add the result to frequency
    %         freq = [freq ff];
    %         if floor(t_closed(dd)/3600) > hour
    %             new_hour(end) = 1;
    %             hour = hour + 1;
    %         end
    %         new_hour = [new_hour 0];
    %         % Reset ff
    %         ff = 1;
    %     end
    % end

    nh = 0;
    % Alternative in case of new hour ==> split the sequence!
    while dd <= length(dir)
        if dir(dd) == dir(dd-1) && floor(t_closed(dd)/3600) == hour
            % We are still in the same direction!
            % Remove direction.
            dir(dd) = [];
            ff = ff + 1;
            t_closed(dd) = [];            
        elseif floor(t_closed(dd)/3600) > hour
            % Change of hour!
            hour = hour + 1;
            freq = [freq ff];
            nh(length(freq)) = 1;

            ff = 0;

        else
            % Change of types!
            % Advance to next spot.
            dd = dd + 1;
            % Add the result to frequency
            freq = [freq ff];
            % Reset ff
            ff = 1;
        end
    end
    % Make sure that freq and new hour are equally long!
    new_hour = zeros(length(freq),1);
    new_hour(1:length(nh)) = nh;

    text = ['\n\t Frequency pattern: '];
    % for ff = 1:length(freq)
    %     if new_hour(ff)
    %         text = [text int2str(freq(ff)) ' */ '];
    %     else
    %         text = [text int2str(freq(ff)) ' / '];
    %     end
    % end

    % Alternative!!!

    text = ['\n\t Frequency pattern: '];
    for ff = 1:length(freq)
        if new_hour(ff)
            text = [text int2str(freq(ff)) ' * '];
        elseif freq(ff) == 0
            text = [text(1:end-1) '/ '];
        else
            text = [text int2str(freq(ff)) ' / '];
        end
    end
    text = [text(1:end-2) '\n'];
    fprintf(text);
    if saveStats
        fprintf(fid, text);
    end
end

% Write the statistics!
text = ['\n***Statistics***\n\n' ...
            '\tCancellations: ' int2str(statistics.nrCancelled) '\n' ...
            '\t\tdirection 0: ' int2str(statistics.nrCancelled_dir0) '\n' ...
            '\t\tdirection 1: ' int2str(statistics.nrCancelled_dir1) '\n'];
if settings.deviation.available
    text = [text ...
            '\tDeviations: ' int2str(statistics.nrDeviated) '\n' ...
            '\t\tdirection 0: ' int2str(statistics.nrDeviated_dir0) '\n' ...
            '\t\tdirection 1: ' int2str(statistics.nrDeviated_dir1) '\n'];
else
    text = [text '\tNo deviation possible\n'];
end
text = [text ...
            '\n\tMaximum delay:\t' statistics.maxDelay_HHMMSS '\n' ...
            '\t\tdirection 0: ' (statistics.maxDelay_HHMMSS_0) '\n' ...
            '\t\tdirection 1: ' (statistics.maxDelay_HHMMSS_1) '\n' ...
            '\tTotal delay:\t' statistics.totalDelay_HHMMSS '\n' ...
            '\tTotal EXTRA delay:\t' statistics.totalExtraDelay_HHMMSS '\n' ...
            '\t\tdirection 0: ' (statistics.totalDelay_HHMMSS_0) '\n' ...
            '\t\tdirection 1: ' (statistics.totalDelay_HHMMSS_1) '\n' ...
            '\tAverage delay:\t' statistics.averageDelay_HHMMSS '\n' ...
            '\t\tdirection 0: ' (statistics.averageDelay_HHMMSS_0) '\n' ...
            '\t\tdirection 1: ' (statistics.averageDelay_HHMMSS_1) '\n'];

fprintf(text);
if saveStats
    fprintf(fid,text);
end

if saveStats
    fclose(fid);
end

try
    measures.orderDirection = dir;
    measures.directionFrequency = freq;
end


statistics.CPUtime = settings.FIFOtime;


if w_dir0 ~= w_dir1
    statistics.OFvalue = Objective;
else 
    statistics.OFvalue = sum(extra_delays) + sum(cancelled)* settings.weights.cancel;
end

    
end

function filename = generateFilename()

filename = 'Measures';
c = clock;

for tt = 1:5
    add = int2str(c(tt));
    if length(add) == 1
        add = ['0' add];
    end
    
    if ismember(tt,[1,4])
        filename = [filename '_'];
    end
    
    filename = [filename add];
end

filename = [filename '.txt'];


end

