function [measures, statistics] = deriveMeasures_v2(timetable, cancelled, releasetime, direction, settings, orig_delays, extra_delays, solu)

saveStats = 1;
% 
% freeOrders_0 = solu.changedOrders_0;
% freeOrders_1 = solu.changedOrders_1;

% Already give some of the statistics here!
% Generate a delay table
delays = [];
row = 0;
trains = unique(timetable.train_id);
for tt = 1:length(trains)
    train_ev = timetable(find(timetable.train_id == trains(tt)),:);
  %  deviated(tt) = train_ev.deviated(1);
end

direction = mod(trains,100);
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

delays_nocancel = [];

if ~isempty(delays)
    delays = struct2table(delays);
    statistics.delays = delays;

    delays_nocancel = orig_delays(find(~cancelled));
    extra_delays_nocancel = extra_delays(find(~cancelled));
    trains_nocancel = trains(find(~cancelled));
    direction_nocancel = direction(find(~cancelled));
end

if ~isempty(delays_nocancel)

    max_delay_rows = find(delays_nocancel == max(delays_nocancel));
    max_delay_row = max_delay_rows(1);
    delays_dir02 = delays_nocancel(direction == 2);
    max_delay_rows_dir02 = find(delays_dir02 == max(delays_dir02));
    max_delay_row_dir02 = max_delay_rows_dir02(1);
    delays_dir03 = delays_nocancel(direction == 3);
    max_delay_rows_dir03 = find(delays_dir03 == max(delays_dir03));
    max_delay_row_dir03 = max_delay_rows_dir03(1);
    delays_dir12 = delays_nocancel(direction == 12);
    max_delay_rows_dir12 = find(delays_dir12 == max(delays_dir12));
    max_delay_row_dir12 = max_delay_rows_dir12(1);
    delays_dir13 = delays_nocancel(direction == 13);
    max_delay_rows_dir13 = find(delays_dir13 == max(delays_dir13));
    max_delay_row_dir13 = max_delay_rows_dir13(1);


    % Max. delay?
    statistics.maxDelay = delays_nocancel(max_delay_row);
    statistics.maxDelay_HHMMSS = timeHHMMSS(delays_nocancel(max_delay_row));
    statistics.maxDelay_train = trains_nocancel(max_delay_rows);
    statistics.maxDelay_02 = delays_dir02(max_delay_row_dir02);
    statistics.maxDelay_HHMMSS_02 = timeHHMMSS(delays_dir02(max_delay_row_dir02));
    statistics.maxDelay_03 = delays_dir03(max_delay_row_dir03);
    statistics.maxDelay_HHMMSS_03 = timeHHMMSS(delays_dir03(max_delay_row_dir03));
    statistics.maxDelay_12 = delays_dir12(max_delay_row_dir12);
    statistics.maxDelay_HHMMSS_12 = timeHHMMSS(delays_dir12(max_delay_row_dir12));
    statistics.maxDelay_13 = delays_dir13(max_delay_row_dir13);
    statistics.maxDelay_HHMMSS_13 = timeHHMMSS(delays_dir13(max_delay_row_dir13));

    % Total delay
    total_delay = sum(delays_nocancel);
    total_extra_delay = sum(extra_delays_nocancel);
    statistics.totalDelay = total_delay;
    statistics.totalExtraDelay = total_extra_delay;
    statistics.totalDelay_HHMMSS = timeHHMMSS(total_delay);
    statistics.totalExtraDelay_HHMMSS = timeHHMMSS(total_extra_delay);
    total_delay_02 = sum(delays_dir02);
    statistics.totalDelay_02 = sum(delays_dir02);
    statistics.totalDelay_HHMMSS_02 = timeHHMMSS(sum(delays_dir02));
    total_delay_03 = sum(delays_dir03);
    statistics.totalDelay_03 = sum(delays_dir03);
    statistics.totalDelay_HHMMSS_03 = timeHHMMSS(sum(delays_dir03));
    total_delay_12 = sum(delays_dir12);
    statistics.totalDelay_12 = sum(delays_dir12);
    statistics.totalDelay_HHMMSS_12 = timeHHMMSS(sum(delays_dir12));
    total_delay_13 = sum(delays_dir13);
    statistics.totalDelay_13 = sum(delays_dir13);
    statistics.totalDelay_HHMMSS_13 = timeHHMMSS(sum(delays_dir13));

    % Average delay
    nr_delay = sum(delays_nocancel > 0);
    av_delay = total_delay / nr_delay;
    statistics.averageDelay = av_delay;
    statistics.averageDelay_HHMMSS = timeHHMMSS(av_delay);
    nr_delay_02 = sum(delays_dir02 > 0);
    av_delay_02 = total_delay_02 / nr_delay_02;
    statistics.averageDelay_02 = av_delay_02;
    statistics.averageDelay_HHMMSS_02 = timeHHMMSS(av_delay_02);
    nr_delay_03 = sum(delays_dir03 > 0);
    av_delay_03 = total_delay_03 / nr_delay_03;
    statistics.averageDelay_03 = av_delay_03;
    statistics.averageDelay_HHMMSS_03 = timeHHMMSS(av_delay_03);
    nr_delay_12 = sum(delays_dir12 > 0);
    av_delay_12 = total_delay_12 / nr_delay_12;
    statistics.averageDelay_12 = av_delay_12;
    statistics.averageDelay_HHMMSS_12 = timeHHMMSS(av_delay_12);
    nr_delay_13 = sum(delays_dir13 > 0);
    av_delay_13 = total_delay_13 / nr_delay_13;
    statistics.averageDelay_13 = av_delay_13;
    statistics.averageDelay_HHMMSS_13 = timeHHMMSS(av_delay_13);
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

%measures.deviatedtrains = trains(find(deviated));
%statistics.deviatedTrains = trains(find(deviated));
%statistics.nrDeviated = sum(deviated);

% Which trains not to consider? The ones after the disruption!
remove = releasetime >= settings.disruption.duration * 3600;
trains(remove) = [];
%t_closed(remove) = [];
cancelled(remove) = [];
direction(remove) = [];

fprintf('***Measures***\n\n');

if saveStats
    fprintf(fid,'***Measures***\n\n');
end

%% First: cancellations
cancelled_dir02 = cancelled(direction == 2);
cancelled_dir03 = cancelled(direction == 3);
cancelled_dir12 = cancelled(direction == 12);
cancelled_dir13 = cancelled(direction == 13);
text02 = [int2str(sum(cancelled)) ' trains have to be cancelled\n \t - direction 02: '...
              int2str(sum(cancelled_dir02))];
text03 = [int2str(sum(cancelled)) ' trains have to be cancelled\n \t - direction 03: '...
              int2str(sum(cancelled_dir03))];
text12 = ['\n\t - direction 12: ' int2str(sum(cancelled_dir12))];
text13 = ['\n\t - direction 13: ' int2str(sum(cancelled_dir13))];

statistics.nrCancelled_dir02 = sum(cancelled_dir02);
statistics.nrCancelled_dir03 = sum(cancelled_dir03);
statistics.nrCancelled_dir12 = sum(cancelled_dir12);
statistics.nrCancelled_dir13 = sum(cancelled_dir13);

text_dir02 = [];
text_dir03 = [];
text_dir12 = [];
text_dir13 = [];
for cc = 1:length(cancelled)
    if cancelled(cc) == 1
        rows = timetable(find(timetable.train_id == trains(cc)),:);
        type = rows.train_type{1};
        time = rows.arrival(1);
        HHMMSS = [int2str(floor(time/3600)) ':' int2str(floor(mod(time,3600))/60) ':' int2str(mod(time,60))];
        text = ['\n\t\t Train ' type int2str(trains(cc)) ' of ' HHMMSS];
        
        if direction(cc) == 2
            text_dir02 = [text_dir02 text];
        end
        if direction(cc) == 3
            text_dir03 = [text_dir03 text];
        end
        if direction(cc) == 12
            text_dir12 = [text_dir12 text];
        end
        if direction(cc) == 13
            text_dir13 = [text_dir13 text];
        end
    end
end
fprintf([text02 text_dir02 text03 text_dir03 text12 text_dir12 text13 text_dir13]);
if saveStats
    fprintf(fid, [text02 text_dir02 text03 text_dir03 text12 text_dir12 text13 text_dir13]);
end
%% If done: deviations
if settings.deviation.available
    deviated_dir0 = deviated(ismember(direction, [2, 3]));
    deviated_dir1 = deviated(ismember(direction, [12, 13]));
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

            if direction(cc) == 12 || direction(cc) == 13
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
% text = [];
% if settings.constraints.fixEntranceOrder
%     text = '\nNo order switching allowed\n';
% elseif (isempty(freeOrders_1) && isempty(freeOrders_0))
%     text = '\nOrder switching allowed, but not done\n';
% else
%     text = '\nChanged orders:\n';
%     if ~isempty(freeOrders_0)
%         text = [text '\t- direction 0: '];
%         for oo = 1:size(freeOrders_0,1)
%             text = [text '(' int2str(freeOrders_0(oo,1)) ', ' ...
%                         int2str(freeOrders_0(oo,2)) '); '];
%         end
%         text = [text '\n'];
%     end
%     if ~isempty(freeOrders_1)
%         text = [text '\t- direction 1: '];
%         for oo = 1:size(freeOrders_1,1)
%             text = [text '(' int2str(freeOrders_1(oo,1)) ', ' ...
%                         int2str(freeOrders_1(oo,2)) '); '];
%         end
%         text = [text '\n'];
%     end
% end
% fprintf(text);
% if saveStats
%     fprintf(fid, text);
% end

% %% Next: for the remaining trains, what order?
% remove = find(cancelled == 1);
% trains(remove) = [];
% t_closed(remove) = [];
% direction(remove) = [];
% 
% % Sort the remaining ones according to their time of entry on the closed
% % segment.
% [t_closed, I] = sortrows(t_closed);
% trains_ordered = trains(I);
% direction_ordered = direction(I);
% 
% if ~isempty(direction_ordered)
%     text = ['\n\n Resulting pattern of trains: '];
%     text_1 = [];    % Train types
%     text_2 = ['\n\t\t\t'];    % Orders
%     for tt = 1:length(I)
%         text_1 = [text_1 int2str(direction_ordered(tt)) ' '];
%         rows = timetable(find(timetable.train_id == trains_ordered(tt)),:);
%         type = rows.train_type{1};
%         text_2 = [text_2 ' ' type];
%     end
%     fprintf([text text_1 text_2]);
%     if saveStats
%         fprintf(fid, [text text_1 text_2]);
%     end
% 
%     % Frequency pattern?
%     text = ['\n\n First direction: ' int2str(direction_ordered(1))];
%     fprintf(text);
%     if saveStats
%         fprintf(fid, text);
%     end
% 
%     dir = direction_ordered;
%     freq = [];
%     dd = 2;
%     ff = 1;
%     hour = 0;
%     new_hour = [];
%     % while dd <= length(dir)
%     %     if dir(dd) == dir(dd-1)
%     %         % We are still in the same direction!
%     %         % Remove direction.
%     %         dir(dd) = [];
%     %         ff = ff + 1;
%     %         if floor(t_closed(dd)/3600) > hour
%     %             % We go to another hour!
%     %             hour = hour + 1;
%     %             new_hour(length(freq)) = 1;
%     %         end
%     %         t_closed(dd) = [];            
%     %     else
%     %         % Change of types!
%     %         % Advance to next spot.
%     %         dd = dd + 1;
%     %         % Add the result to frequency
%     %         freq = [freq ff];
%     %         if floor(t_closed(dd)/3600) > hour
%     %             new_hour(end) = 1;
%     %             hour = hour + 1;
%     %         end
%     %         new_hour = [new_hour 0];
%     %         % Reset ff
%     %         ff = 1;
%     %     end
%     % end
% 
%     nh = 0;
%     % Alternative in case of new hour ==> split the sequence!
%     while dd <= length(dir)
%         if dir(dd) == dir(dd-1) && floor(t_closed(dd)/3600) == hour
%             % We are still in the same direction!
%             % Remove direction.
%             dir(dd) = [];
%             ff = ff + 1;
%             t_closed(dd) = [];            
%         elseif floor(t_closed(dd)/3600) > hour
%             % Change of hour!
%             hour = hour + 1;
%             freq = [freq ff];
%             nh(length(freq)) = 1;
% 
%             ff = 0;
% 
%         else
%             % Change of types!
%             % Advance to next spot.
%             dd = dd + 1;
%             % Add the result to frequency
%             freq = [freq ff];
%             % Reset ff
%             ff = 1;
%         end
%     end
%     % Make sure that freq and new hour are equally long!
%     new_hour = zeros(length(freq),1);
%     new_hour(1:length(nh)) = nh;
% 
%     text = ['\n\t Frequency pattern: '];
%     % for ff = 1:length(freq)
%     %     if new_hour(ff)
%     %         text = [text int2str(freq(ff)) ' */ '];
%     %     else
%     %         text = [text int2str(freq(ff)) ' / '];
%     %     end
%     % end
% 
%     % Alternative!!!
% 
%     text = ['\n\t Frequency pattern: '];
%     for ff = 1:length(freq)
%         if new_hour(ff)
%             text = [text int2str(freq(ff)) ' * '];
%         elseif freq(ff) == 0
%             text = [text(1:end-1) '/ '];
%         else
%             text = [text int2str(freq(ff)) ' / '];
%         end
%     end
%     text = [text(1:end-2) '\n'];
%     fprintf(text);
%     if saveStats
%         fprintf(fid, text);
%     end
% end
% 
% 
% 
% 
% statistics.CPUtime = solu.solvertime;
% try   % For CPLEX
%     statistics.OFvalue = solu.solveroutput.fval;
% catch % For GUROBI
%     statistics.OFvalue = solu.solveroutput.result.objval;
% end
    
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

