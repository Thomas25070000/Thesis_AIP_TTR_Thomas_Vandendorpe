function [] = boxplotDelays(input, plottitle, splitDirections)

if nargin < 3
    splitDirection = 1;
end


points = [];
counter = 0;
for tt = 1:size(input,2)
    TT = input(tt).timetables;
    delays = input(tt).statistics.delays;
    
    ll = input(tt).label;
    if iscell(ll)
        ll = ll{:};
    end

    for dd = 1:size(delays,1)
        
        counter = counter + 1;

        train = delays.train_id(dd);
        dir(counter) = mod(train,10);
        value = delays.orig_delay(dd);
        
        if splitDirections
            points(counter).label = sprintf([ll ' (' int2str(dir(counter)) ')']);
        else
            points(counter).label = sprintf([ll]);
        end
            
        % Departure time?
        train_ev = TT(find(TT.train_id == train),:);
        deptime = train_ev.arrival(1);
        type = train_ev.train_type{1};

        points(counter).X = deptime;
        points(counter).Y = value;
        if dir
            points(counter).symbol = 's';
            points(counter).color = 'blue';
        else
            points(counter).symbol = 'o';
            points(counter).color = 'black';
        end
    end
end

% Plot the values
figure;
hold on

boxplot([points.Y], {points.label})
title(plottitle);
ylabel('Delays [s]')

hold off




end