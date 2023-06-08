function [timetable, blocksections, runningtimes] = slimTimetable(timetable,blocksections,runningtimes)

% Remove the blocksections A and B!
AandB = find(ismember(blocksections.type,{'A','B'}));
blocksections(AandB,:) = [];
timetable(ismember(timetable.blocksection,AandB),:) = [];

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





runningtimes.L1.disrupted(AandB) = [];
runningtimes.L1.regular(AandB) = [];

runningtimes.L0.disrupted(AandB) = [];
runningtimes.L0.regular(AandB) = [];

runningtimes.IC1.disrupted(AandB) = [];
runningtimes.IC1.regular(AandB) = [];

runningtimes.IC0.disrupted(AandB) = [];
runningtimes.IC0.regular(AandB) = [];

end