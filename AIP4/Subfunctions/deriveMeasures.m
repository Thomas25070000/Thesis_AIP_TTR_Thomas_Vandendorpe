function measures = deriveMeasures(timetable, t_closed, cancelled, releasetime, direction, settings)

% Already give some of the statistics here!


folder = [pwd '\Output'];

filename = generateFilename();
folderfile = [folder '\' filename];

fid=fopen(folderfile,'wt');
if fid <= 0
    % Problem with opening the file! Mostly because there is no folder
    % present.
    mkdir('Output');
    fid=fopen(folderfile,'wt');
end

trains = unique(timetable.train_id);

measures.cancelledtrains = trains(find(cancelled));

% Which trains not to consider? The ones after the disruption!
remove = releasetime >= settings.disruption.duration * 3600;
trains(remove) = [];
t_closed(remove) = [];
cancelled(remove) = [];
direction(remove) = [];

%% First: cancellations
cancelled_dir0 = cancelled(direction == 0);
cancelled_dir1 = cancelled(direction == 1);
text0 = [int2str(sum(cancelled)) ' trains have to be cancelled\n \t - direction 0: '...
              int2str(sum(cancelled_dir0))];
text1 = ['\n\t - direction 1: ' int2str(sum(cancelled_dir1))];

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
fprintf(fid, [text0 text_dir0 text1 text_dir1]);


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
fprintf(fid, [text text_1 text_2]);

% Frequency pattern?
text = ['\n\n First direction: ' int2str(direction_ordered(1))];
fprintf(text);
fprintf(fid, text);

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
fprintf(text);
fprintf(fid, text);

fclose(fid);

measures.orderDirection = dir;
measures.directionFrequency = freq;
    
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