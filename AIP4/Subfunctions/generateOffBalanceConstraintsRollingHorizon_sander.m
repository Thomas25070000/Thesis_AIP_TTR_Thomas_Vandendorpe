function balanceConstraints = generateOffBalanceConstraintsRollingHorizon(trains, releasetime, x, settings)

delta = settings.disruption.offbalance;

% In which direction do they go?
direction = mod(trains,10);

dir0 = find(direction == 0);
dir1 = find(direction == 1);

release0 = releasetime(dir0);
release1 = releasetime(dir1);

% Frequency of trains?
freq = length(trains) / (2*(settings.disruption.duration + settings.disruption.nr_hours_after));

if min(release1) <= min(release0)
    C = [];
    for tt = 1:length(dir1)-freq+1
        % Which run in direction 0 in this interval?
        time = release1(tt);
        t0 = dir0(find(release0 >= time & release0 < (time + 3600 - 1)));
        t1 = dir1(find(release1 >= time & release1 < (time + 3600 - 1)));

        LB = [sum(x((t0))) - sum(x((t1))) >= -delta];
        UB = [sum(x((t0))) - sum(x((t1))) <= delta];

        C = [C LB UB];
    end
else
    C = [];
    for tt = 1:length(dir0)-freq+1
        % Which run in direction 0 in this interval?
        time = release0(tt);
        t0 = dir0(find(release0 >= time & release0 < (time + 3600 - 1)));
        t1 = dir1(find(release1 >= time & release1 < (time + 3600 - 1)));

        LB = [sum(x((t0))) - sum(x((t1))) >= -delta];
        UB = [sum(x((t0))) - sum(x((t1))) <= delta];

        C = [C LB UB];
    end
end

balanceConstraints = C;



end