function consecutiveConstraints = generateConsecutiveConstraints(timetable, trains, settings, q, x)

limit = settings.TT.orderlimit;

Ntrains = length(trains);

% Get release times again.
releasetime = zeros(Ntrains,1);
for tt = 1:Ntrains
    rows_train = find(timetable.train_id == trains(tt));
    releasetime(tt) = timetable.arrival(rows_train(1));
end

% The ones after the disruption are not considered.
remove = releasetime >= settings.disruption.duration * 3600;
trains(remove) = [];
releasetime(remove) = [];

% Sort trains per direction, and according to increasing release time.
dir0 = find(mod(trains,10) == 0);
dir1 = find(mod(trains,10) == 1);

release0 = releasetime(dir0);
release1 = releasetime(dir1);

[~, order0] = sort(release0);
[~, order1] = sort(release1);

dir0_ordered = dir0(order0);
dir1_ordered = dir1(order1);    

% For both directions, generate the constraints for each pair.
Constraints = [];
for ii = 2:length(dir0_ordered)
    % Which trains / q values do we need?
    first = dir0_ordered(ii-1);
    second = dir0_ordered(ii);
    
    % Compare with all trains of the direction, except the same one.
%     first_comp = dir0(dir0 ~= first);
%     second_comp = dir0(dir0 ~= second);
    first_comp = dir1;
    second_comp = dir1;
    
    label = ['nr_trains_between_' int2str(trains(first)) '_' int2str(trains(second))];
    
    Q1 = 0;
    for tt = 1:length(first_comp)
        if first_comp(tt) < first
            Q1 = Q1 + q(first_comp(tt),first);
        else
            Q1 = Q1 + 1 - q(first,first_comp(tt));
        end
    end
    Q2 = 0;
    for tt = 1:length(second_comp)
        if second_comp(tt) < second
            Q2 = Q2 + q(second_comp(tt),second);
        else
            Q2 = Q2 + 1 - q(second,second_comp(tt));
        end
    end
    
    C = [Q2 - Q1 <= limit]: label;
    Constraints = [Constraints, C];
end

for ii = 2:length(dir1_ordered)
    % Which trains / q values do we need?
    first = dir1_ordered(ii-1);
    second = dir1_ordered(ii);
    
    % Compare with all trains of the direction, except the same one.
%     first_comp = dir1(dir1 ~= first);
%     second_comp = dir1(dir1 ~= second);
    first_comp = dir0;
    second_comp = dir0;
    
    label = ['nr_trains_between_' int2str(trains(first)) '_' int2str(trains(second))];
    
    Q1 = 0;
    for tt = 1:length(first_comp)
        if first_comp(tt) < first
            Q1 = Q1 + q(first_comp(tt),first);
        elseif first_comp(tt) > first
            Q1 = Q1 + 1 - q(first,first_comp(tt));
        end
    end
    Q2 = 0;
    for tt = 1:length(second_comp)
        if second_comp(tt) < second
            Q2 = Q2 + q(second_comp(tt),second);
        elseif second_comp(tt) > second
            Q2 = Q2 + 1 - q(second,second_comp(tt));
        end
    end
    
    C = [Q2 - Q1 <= limit]: label;
    Constraints = [Constraints, C];
end
    
consecutiveConstraints = Constraints;
    
    
    
    
    
end