function [measures, statistics] = deriveMeasures_case3(timetable, traininfo, settings, solu);

saveStats = settings.saveStats;

% Already give some of the statistics here!
% Generate a delay table
delays = [];
row = 0;
trains = [traininfo.id];
direction = [traininfo.dir];

orig_delays = [traininfo.delay];
extra_delays = [traininfo.extra_delay];
cancelled = [traininfo.cancelled];
releasetime = [traininfo.entry];
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

    delays_nocancel = orig_delays(find(cancelled == 0));
    extra_delays_nocancel = extra_delays(find(cancelled == 0));
    trains_nocancel = trains(find(cancelled == 0));
    direction_nocancel = direction(find(cancelled == 0));
else
    delays_nocancel = [];
end

if ~isempty(delays_nocancel)

    max_delay_rows = find(delays_nocancel == max(delays_nocancel));
    max_delay_row = max_delay_rows(1);
    delays_dir02 = delays_nocancel(find(direction_nocancel == 12));
    max_delay_rows_dir02 = delays_dir02 == max(delays_dir02);
    delays_dir03 = delays_nocancel(find(direction_nocancel == 13));
    max_delay_rows_dir03 = delays_dir03 == max(delays_dir03);
    try
        max_delay_row_dir02 = max_delay_rows_dir02(1);
        statistics.maxDelay_02 = delays_dir02(max_delay_row_dir02);
        statistics.maxDelay_HHMMSS_02 = timeHHMMSS(delays_dir02(max_delay_row_dir02));
    catch
        max_delay_row_dir02 = 0;
        statistics.maxDelay_02 = 0;
        statistics.maxDelay_HHMMSS_02 = 0;
    end
    try
        max_delay_row_dir03 = max_delay_rows_dir03(1);
        statistics.maxDelay_03 = delays_dir03(max_delay_row_dir03);
        statistics.maxDelay_HHMMSS_03 = timeHHMMSS(delays_dir03(max_delay_row_dir03));
    catch
        max_delay_row_dir03 = 0;
        statistics.maxDelay_03 = 0;
        statistics.maxDelay_HHMMSS_03 = 0;
    end
    delays_dir12 = delays_nocancel(find(direction_nocancel == 12));
    max_delay_rows_dir12 = find(delays_dir12 == max(delays_dir12));
    try
        max_delay_row_dir12 = max_delay_rows_dir12(1);
        statistics.maxDelay_12 = delays_dir12(max_delay_row_dir12);
        statistics.maxDelay_HHMMSS_12 = timeHHMMSS(delays_dir12(max_delay_row_dir12));
    catch
        max_delay_row_dir12 = 0;
        statistics.maxDelay_12 = 0;
        statistics.maxDelay_HHMMSS_12 = 0;
    end
    delays_dir13 = delays_nocancel(find(direction_nocancel == 13));
    max_delay_rows_dir13 = find(delays_dir13 == max(delays_dir13));
    try
        max_delay_row_dir13 = max_delay_rows_dir13(1);
        statistics.maxDelay_13 = delays_dir13(max_delay_row_dir13);
        statistics.maxDelay_HHMMSS_13 = timeHHMMSS(delays_dir13(max_delay_row_dir13));
    catch
        max_delay_row_dir13 = 0;
        statistics.maxDelay_13 = 0;
        statistics.maxDelay_HHMMSS_13 = 0;
    end

    % Max. delay?
    statistics.maxDelay = delays_nocancel(max_delay_row);
    statistics.maxDelay_HHMMSS = timeHHMMSS(delays_nocancel(max_delay_row));
    statistics.maxDelay_train = trains_nocancel(max_delay_rows);
%     statistics.maxDelay_0 = delays_dir0(max_delay_row_dir0);
%     statistics.maxDelay_HHMMSS_0 = timeHHMMSS(delays_dir0(max_delay_row_dir0));
%     statistics.maxDelay_1 = delays_dir1(max_delay_row_dir1);
%     statistics.maxDelay_HHMMSS_1 = timeHHMMSS(delays_dir1(max_delay_row_dir1));

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
    statistics.maxDelay_02 = 0;
    statistics.maxDelay_HHMMSS_02 = 0;
    statistics.maxDelay_03 = 0;
    statistics.maxDelay_HHMMSS_02 = 0;
    statistics.maxDelay_1 = 0;
    statistics.maxDelay_HHMMSS_1 = 0;
    statistics.maxDelay_12 = 0;
    statistics.maxDelay_HHMMSS_13 = 0;
    statistics.maxDelay_train = 0;
    statistics.totalDelay = 0;
    statistics.totalDelay_HHMMSS = 0;
    statistics.totalExtraDelay = 0;
    statistics.totalExtraDelay_HHMMSS = 0;
    statistics.totalDelay_02 = 0;
    statistics.totalDelay_HHMMSS_02 = 0;
    statistics.totalDelay_03 = 0;
    statistics.totalDelay_HHMMSS_03 = 0;
    statistics.totalDelay_12 = 0;
    statistics.totalDelay_HHMMSS_12 = 0;
    statistics.totalDelay_13 = 0;
    statistics.totalDelay_HHMMSS_13 = 0;
    statistics.averageDelay = 0;
    statistics.averageDelay_HHMMSS = 0;
    statistics.averageDelay_02 = 0;
    statistics.averageDelay_HHMMSS_02 = 0;
    statistics.averageDelay_03 = 0;
    statistics.averageDelay_HHMMSS_03 = 0;
    statistics.averageDelay_12 = 0;
    statistics.averageDelay_HHMMSS_12 = 0;
    statistics.averageDelay_13 = 0;
    statistics.averageDelay_HHMMSS_13 = 0;
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

measures.cancelledtrains = trains(find(cancelled));

statistics.cancelledTrains = trains(find(cancelled));
statistics.nrCancelled = sum(cancelled);


% Which trains not to consider? The ones after the disruption!
remove = releasetime >= settings.disruption.duration * 3600;
trains(remove) = [];
cancelled(remove) = [];
direction(remove) = [];

fprintf('***Measures***\n\n');

if saveStats
    fprintf(fid,'***Measures***\n\n');
end

%% First: cancellations
cancelled_dir02 = cancelled(direction == 02);
cancelled_dir03 = cancelled(direction == 03);
cancelled_dir12 = cancelled(direction == 12);
cancelled_dir13 = cancelled(direction == 13);
% text0 = [int2str(sum(cancelled)) ' trains have to be cancelled\n \t - direction 0: '...
%               int2str(sum(cancelled_dir0))];
% text1 = ['\n\t - direction 1: ' int2str(sum(cancelled_dir1))];

statistics.nrCancelled_dir02 = sum(cancelled_dir02);
statistics.nrCancelled_dir03 = sum(cancelled_dir03);
statistics.nrCancelled_dir12 = sum(cancelled_dir12);
statistics.nrCancelled_dir13 = sum(cancelled_dir13);
% 
% text_dir0 = [];
% text_dir1 = [];
% for cc = 1:length(cancelled)
%     if cancelled(cc) == 1
%         type = traininfo(cc).type;
%         time = traininfo(cc).entry;
%         HHMMSS = [int2str(floor(time/3600)) ':' int2str(floor(mod(time,3600))/60) ':' int2str(mod(time,60))];
%         text = ['\n\t\t Train ' type int2str(trains(cc)) ' of ' HHMMSS];
%         
%         if direction(cc) == 1
%             text_dir1 = [text_dir1 text];
%         else
%             text_dir0 = [text_dir0 text];
%         end
%     end
% end
% fprintf([text0 text_dir0 text1 text_dir1]);
% if saveStats
%     fprintf(fid, [text0 text_dir0 text1 text_dir1]);
% end



% %% Order switching?
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

%% Next: for the remaining trains, on which track?
% remove = find(cancelled == 1);
% trains(remove) = [];
% direction(remove) = [];
% 
% % For each track
% text = ['\nScheduled on tracks\n'];
% 
% Ntracks = size(settings.tracks,1);
% measures_tracks = settings.tracks;
% measures_tracks.dir0(1) = 0;
% measures_tracks.dir1(1) = 0;
% for tr = 1:Ntracks
%     whichtrains = find([traininfo.newtrack] == tr);
%     track_dir0 = whichtrains(find([traininfo(whichtrains).dir] == 0));
%     track_dir1 = whichtrains(find([traininfo(whichtrains).dir] == 1));
%     speed = settings.tracks.vdis(tr);
%     text = [text '\t- track ' int2str(tr) ' (speed = ' int2str(speed) 'km/h)\n'];
%     text = [text '\t\tdir 0: ' int2str(length(track_dir0)) '\n'];
%     text = [text '\t\tdir 1: ' int2str(length(track_dir1)) '\n'];
%     
%     measures_tracks.dir0(tr) = length(track_dir0);
%     measures_tracks.dir1(tr) = length(track_dir1);
% end
% measures.tracks = measures_tracks;
% 
% fprintf(text);
% if saveStats
%     fprintf(fid,text);
% end




% %% Write the statistics!
% text = ['\n***Statistics***\n\n' ...
%             '\tCancellations: ' int2str(statistics.nrCancelled) '\n' ...
%             '\t\tdirection 0: ' int2str(statistics.nrCancelled_dir0) '\n' ...
%             '\t\tdirection 1: ' int2str(statistics.nrCancelled_dir1) '\n'];
% if settings.deviation.available
%     text = [text ...
%             '\tDeviations: ' int2str(statistics.nrDeviated) '\n' ...
%             '\t\tdirection 0: ' int2str(statistics.nrDeviated_dir0) '\n' ...
%             '\t\tdirection 1: ' int2str(statistics.nrDeviated_dir1) '\n'];
% else
%     text = [text '\tNo deviation possible\n'];
% end
% text = [text ...
%             '\n\tMaximum delay:\t' statistics.maxDelay_HHMMSS '\n' ...
%             '\t\tdirection 0: ' (statistics.maxDelay_HHMMSS_0) '\n' ...
%             '\t\tdirection 1: ' (statistics.maxDelay_HHMMSS_1) '\n' ...
%             '\tTotal delay:\t' statistics.totalDelay_HHMMSS '\n' ...
%             '\tTotal EXTRA delay:\t' statistics.totalExtraDelay_HHMMSS '\n' ...
%             '\t\tdirection 0: ' (statistics.totalDelay_HHMMSS_0) '\n' ...
%             '\t\tdirection 1: ' (statistics.totalDelay_HHMMSS_1) '\n' ...
%             '\tAverage delay:\t' statistics.averageDelay_HHMMSS '\n' ...
%             '\t\tdirection 0: ' (statistics.averageDelay_HHMMSS_0) '\n' ...
%             '\t\tdirection 1: ' (statistics.averageDelay_HHMMSS_1) '\n'];
% 
% fprintf(text);
% if saveStats
%     fprintf(fid,text);
% end
% 
% if saveStats
%     fclose(fid);
% end



statistics.CPUtime = solu.solvertime;
try   % For CPLEX
    statistics.OFvalue = solu.solveroutput.fval;
catch % For GUROBI
    statistics.OFvalue = solu.solveroutput.result.objval;
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

