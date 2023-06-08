function timetable = delayTrain(timetable, train_id, delay, complete)

% Do we want to get the complete timetable back as return?
% complete = 1 (yes) or = 0 (no)


if isempty(train_id)
    train_id = unique(timetable.train_id);
end

rows = find(ismember(timetable.train_id, train_id));

% Delay it
timetable.arrival(rows) = timetable.arrival(rows) + delay;
timetable.departure(rows) = timetable.departure(rows) + delay;
timetable.start(rows) = timetable.start(rows) + delay;
timetable.finish(rows) = timetable.finish(rows) + delay;

if ~complete
    timetable = timetable(rows,:);
end


end