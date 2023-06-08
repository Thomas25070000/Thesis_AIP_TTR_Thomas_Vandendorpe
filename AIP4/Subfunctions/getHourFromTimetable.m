function hour_tt = getHourFromTimetable(timetable,settings,hour)

hour_tt = [];
maxHour = settings.disruption.duration + settings.disruption.nr_hours_after;

if nargin < 3
    % Hour has not been given, just return the most busy one!
    hh = floor([unique(timetable.train_id)]/1000);
    N = hist(hh,max(hh + 1));
    HHreturn = find(N == max(N));
    HHreturn = HHreturn(1);
    rows = find(floor([timetable.train_id]/1000) == HHreturn - 1);
    hour_tt = timetable(rows,:);
elseif hour <= maxHour
    % Get the rows
    rows = find(floor([timetable.train_id]/1000) == hour - 1);
    hour_tt = timetable(rows,:);
else
    disp('No valid hour!');
    pause;
end










end