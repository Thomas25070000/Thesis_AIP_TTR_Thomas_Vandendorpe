function [timetable, blocksections, runningtimes] = slimTimetable(timetable,blocksections,runningtimes)

% Remove the blocksections A and B!
AandBandC = find(ismember(blocksections.type,{'A','B','C'}));
blocksections(AandBandC,:) = [];
timetable(ismember(timetable.blocksection,AandBandC),:) = [];

% Replace the indices
replace = unique(blocksections.id);
blocksections.id(:) = 1:length(replace);

for ee = 1:size(timetable,1)
    timetable.blocksection(ee) = find(timetable.blocksection(ee)==replace);
end

timetable.orig_id(1) = 0;
for ee = 1:size(timetable,1)
    timetable.orig_id(ee) = timetable.event_id(ee);
    timetable.event_id(ee) = ee;
end





runningtimes.L12.disrupted(AandBandC) = [];
runningtimes.L12.regular(AandBandC) = [];

runningtimes.L13.disrupted(AandBandC) = [];
runningtimes.L13.regular(AandBandC) = [];

runningtimes.L02.disrupted(AandBandC) = [];
runningtimes.L02.regular(AandBandC) = [];

runningtimes.L03.disrupted(AandBandC) = [];
runningtimes.L03.regular(AandBandC) = [];

runningtimes.IC12.disrupted(AandBandC) = [];
runningtimes.IC12.regular(AandBandC) = [];

runningtimes.IC13.disrupted(AandBandC) = [];
runningtimes.IC13.regular(AandBandC) = [];

runningtimes.IC02.disrupted(AandBandC) = [];
runningtimes.IC02.regular(AandBandC) = [];

runningtimes.IC03.disrupted(AandBandC) = [];
runningtimes.IC03.regular(AandBandC) = [];

end