function consecutiveConstraints = generateConsecutiveConstraints_v3(timetable, trains, settings, q, x, dev)

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
dir12 = find(mod(trains,10) == 2);
dir13 = find(mod(trains,10) == 3);

release0 = releasetime(dir0);
release12 = releasetime(dir12);
release13 = releasetime(dir13);


[release0, order0] = sort(release0);
[release12, order12] = sort(release12);
[release13, order13] = sort(release13);


dir0_ordered = dir0(order0);
dir12_ordered = dir12(order12);    
dir13_ordered = dir13(order13);    


% For both directions, generate the constraints for each pair.
Constraints = [];
for ii = 1:length(dir0_ordered)-1
    for jj = (ii+1):length(dir0_ordered)
        % Which trains / q values do we need?
        first = dir0_ordered(ii);
        second = dir0_ordered(jj);

        % Compare with all trains of the direction, except the same one.
    %     first_comp = dir0(dir0 ~= first);
    %     second_comp = dir0(dir0 ~= second);
        first_comp_12 = dir12;
        second_comp_12 = dir12;
        first_comp_13 = dir13;
        second_comp_13 = dir13;

        label = ['nr_trains_between_' int2str(trains(first)) '_' int2str(trains(second))];

        % Are there trains inbetween?
        inbetween = dir0_ordered((ii+1):(jj-1));
        % If it is empty, set cancel part to 0
        try
            cancelpart = sum(x(inbetween) + dev(inbetween));
        catch
            cancelpart = 0;
        end
        % Derive the "cancel factor"
        cancel_factor = jj - ii - cancelpart;
        
        % Create the sums to count.
        Q1 = 0;
        for tt = 1:length(first_comp_12)
            if first_comp_12(tt) < first
                Q1 = Q1 + q(first_comp_12(tt),first);
            else
                Q1 = Q1 + 1 - q(first,first_comp_12(tt));
            end
        end
        for tt = 1:length(first_comp_13)
            if first_comp(tt) < first
                Q1 = Q1 + q(first_comp_13(tt),first);
            else
                Q1 = Q1 + 1 - q(first,first_comp_13(tt));
            end
        end
        Q2 = 0;
        for tt = 1:length(second_comp_12)
            if second_comp_12(tt) < second
                Q2 = Q2 + q(second_comp_12(tt),second);
            else
                Q2 = Q2 + 1 - q(second,second_comp_12(tt));
            end
        end
        for tt = 1:length(second_comp_13)
            if second_comp_13(tt) < second
                Q2 = Q2 + q(second_comp_13(tt),second);
            else
                Q2 = Q2 + 1 - q(second,second_comp_13(tt));
            end
        end

        C = [Q2 - Q1 - limit * cancel_factor - (length(second_comp_12)+length(second_comp_13)) * (x(first) + dev(first) + x(second) + dev(second) + gamma(first)) <= 0]: label;
        Constraints = [Constraints, C];
    end
    
    % Which ones before the first one?
   	Q = 0;
    first = dir0_ordered(ii);
    first_comp_12 = dir12;
    first_comp_13 = dir13;

    for tt = 1:length(first_comp_12)
        if first_comp_12(tt) < first
            Q = Q + q(first_comp_12(tt),first);
        else
            Q = Q + 1 - q(first,first_comp_12(tt));
        end
    end 
    for tt = 1:length(first_comp_13)
        if first_comp_13(tt) < first
            Q = Q + q(first_comp_13(tt),first);
        else
            Q = Q + 1 - q(first,first_comp_13(tt));
        end
    end 
    % Are there trains inbetween?
    inbetween = dir0_ordered(1:(ii-1));
    % If it is empty, set cancel part to 0
    try
        cancelpart = sum(x(inbetween) + dev(inbetween));
    catch
        cancelpart = 0;
    end
    % Derive the "cancel factor"
    cancel_factor = ii - cancelpart;
    
    firstC = [Q - sum(x(dir12) + dev(dir12) + x(dir13) + dev(dir13)) - limit * cancel_factor - (length(second_comp_12)+length(second_comp_13)) * (x(first) + dev(first)) <= 0]: label;
    Constraints = [Constraints, firstC];
    
end

for ii = 1:length(dir12_ordered)-1
    for jj = (ii+1):length(dir12_ordered)
        % Which trains / q values do we need?
        first = dir12_ordered(ii);
        second = dir12_ordered(jj);

        % Compare with all trains of the direction, except the same one.
    %     first_comp = dir0(dir0 ~= first);
    %     second_comp = dir0(dir0 ~= second);
        first_comp = dir0;
        second_comp = dir0;

        label = ['nr_trains_between_' int2str(trains(first)) '_' int2str(trains(second))];

        % Are there trains inbetween?
        inbetween = dir12_ordered((ii+1):(jj-1));
        % If it is empty, set cancel part to 0
        try
            cancelpart = sum(x(inbetween) + dev(inbetween));
        catch
            cancelpart = 0;
        end
        % Derive the "cancel factor"
        cancel_factor = jj - ii - cancelpart;
        
        % Create the sums to count.
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

        C = [Q2 - Q1 - limit * cancel_factor - length(second_comp) * (x(first) + dev(first) + x(second) + dev(second) + gamma(first)) <= 0]: label;
        Constraints = [Constraints, C];
    end
    
    % Which ones before the first one?
   	Q = 0;
    first = dir12_ordered(ii);
    first_comp = dir0;
    for tt = 1:length(first_comp)
        if first_comp(tt) < first
            Q = Q + q(first_comp(tt),first);
        else
            Q = Q + 1 - q(first,first_comp(tt));
        end
    end 
    % Are there trains inbetween?
    inbetween = dir12_ordered(1:(ii-1));
    % If it is empty, set cancel part to 0
    try
        cancelpart = sum(x(inbetween) + dev(inbetween));
    catch
        cancelpart = 0;
    end
    % Derive the "cancel factor"
    cancel_factor = ii - cancelpart;
    
    firstC = [Q - sum(x(dir0) + dev(dir0)) - limit * cancel_factor - length(second_comp) * (x(first) + dev(first)) <= 0]: label;
    Constraints = [Constraints, firstC];
    
    
end

for ii = 1:length(dir13_ordered)-1
    for jj = (ii+1):length(dir13_ordered)
        % Which trains / q values do we need?
        first = dir13_ordered(ii);
        second = dir13_ordered(jj);

        % Compare with all trains of the direction, except the same one.
    %     first_comp = dir0(dir0 ~= first);
    %     second_comp = dir0(dir0 ~= second);
        first_comp = dir0;
        second_comp = dir0;

        label = ['nr_trains_between_' int2str(trains(first)) '_' int2str(trains(second))];

        % Are there trains inbetween?
        inbetween = dir13_ordered((ii+1):(jj-1));
        % If it is empty, set cancel part to 0
        try
            cancelpart = sum(x(inbetween) + dev(inbetween));
        catch
            cancelpart = 0;
        end
        % Derive the "cancel factor"
        cancel_factor = jj - ii - cancelpart;
        
        % Create the sums to count.
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

        C = [Q2 - Q1 - limit * cancel_factor - length(second_comp) * (x(first) + dev(first) + x(second) + dev(second) + gamma(first)) <= 0]: label;
        Constraints = [Constraints, C];
    end
    
    % Which ones before the first one?
   	Q = 0;
    first = dir13_ordered(ii);
    first_comp = dir0;
    for tt = 1:length(first_comp)
        if first_comp(tt) < first
            Q = Q + q(first_comp(tt),first);
        else
            Q = Q + 1 - q(first,first_comp(tt));
        end
    end 
    % Are there trains inbetween?
    inbetween = dir13_ordered(1:(ii-1));
    % If it is empty, set cancel part to 0
    try
        cancelpart = sum(x(inbetween) + dev(inbetween));
    catch
        cancelpart = 0;
    end
    % Derive the "cancel factor"
    cancel_factor = ii - cancelpart;
    
    firstC = [Q - sum(x(dir0) + dev(dir0)) - limit * cancel_factor - length(second_comp) * (x(first) + dev(first)) <= 0]: label;
    Constraints = [Constraints, firstC];
    
    
end

    
consecutiveConstraints = Constraints;
    
    
    
    
    
end