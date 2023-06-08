function [arrtime, arrtimeHHMMSS, deptime, deptimeHHMMSS] = retrieveArrivalAndDepartureTimes(timetable, firstTime)

trains = unique(timetable.train_id);

for tt = 1:length(trains)
    events = timetable(find(timetable.train_id == trains(tt)),:);
    arr(tt) = events.arrival(1);
    dep(tt) = events.departure(end);
end

arr = arr + firstTime;
dep = dep + firstTime;

dir02 = find(mod(trains,100) == 2);
dir03 = find(mod(trains,100) == 3);
dir12 = find(mod(trains,100) == 12);
dir13 = find(mod(trains,100) == 13);
[~, I12] = sort(arr(dir12));
[~, I13] = sort(arr(dir13));
[~, I02] = sort(arr(dir02));
[~, I03] = sort(arr(dir03));

order = [dir12(I12); dir13(I13); dir02(I02); dir03(I03)];

deptime = dep(order)';
arrtime = arr(order)';

deptimeHHMMSS = {};
for dd = 1:length(deptime)
    deptimeHHMMSS{dd} = timeHHMMSS(deptime(dd));
end
for aa = 1:length(arrtime)
    arrtimeHHMMSS{aa} = timeHHMMSS(arrtime(aa));
end
deptimeHHMMSS = deptimeHHMMSS';
arrtimeHHMMSS = arrtimeHHMMSS';


end