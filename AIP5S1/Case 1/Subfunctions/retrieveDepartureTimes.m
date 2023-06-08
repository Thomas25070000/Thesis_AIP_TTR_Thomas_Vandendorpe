function [deptime, deptimeHHMMSS] = retrieveDepartureTimes(timetable, firstTime)

trains = unique(timetable.train_id);

for tt = 1:length(trains)
    events = timetable(find(timetable.train_id == trains(tt)),:);
    arr(tt) = events.arrival(1);
    dep(tt) = events.departure(end);
end

arr = arr + firstTime;
dep = dep + firstTime;

dir0 = find(mod(trains,10) == 0);
dir1 = find(mod(trains,10) == 1);
[~, I1] = sort(arr(dir1));
[~, I0] = sort(arr(dir0));

order = [dir1(I1); dir0(I0)];

deptime = dep(order)';

deptimeHHMMSS = {};
for dd = 1:length(deptime)
    deptimeHHMMSS{dd} = timeHHMMSS(deptime(dd));
end
deptimeHHMMSS = deptimeHHMMSS';





end