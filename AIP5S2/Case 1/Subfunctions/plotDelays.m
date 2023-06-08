function [] = plotDelays(new_timetable, statistics,settings)

TT = new_timetable;
delays = statistics.delays;

points = [];
for dd = 1:size(delays,1)
    
    train = delays.train_id(dd);
    dir(dd) = mod(train,10);
    value = delays.orig_delay(dd);
    
    % Departure time?
    train_ev = TT(find(TT.train_id == train),:);
    deptime = train_ev.arrival(1);
    type = train_ev.train_type{1};
    
    points(dd).X = deptime;
    points(dd).Y = value;
    if dir
        points(dd).symbol = 's';
        points(dd).color = 'blue';
    else
        points(dd).symbol = 'o';
        points(dd).color = 'black';
    end
    
%     switch type
%         case 'IC'
%             points(dd).color = 'black';
%         case 'R'
%             points(dd).color = 'blue';
%             points(dd).color = 'black';
%     end
end


% Plot the values
figure;
hold on

% First, direction 0
X = [];
Y = [];
for pp = 1:size(points,2)
    if dir(pp) == 0
        X = [X points(pp).X];
        Y = [Y points(pp).Y];
    end
end
dot = scatter(X,Y);
dot.Marker = 'o';
dot.MarkerFaceColor = 'black';
dot.MarkerEdgeColor = 'none'; 


X = [];
Y = [];
for pp = 1:size(points,2)
    if dir(pp) == 1
        X = [X points(pp).X];
        Y = [Y points(pp).Y];
    end
end
dot = scatter(X,Y);
dot.Marker = 's';
dot.MarkerFaceColor = 'blue';
dot.MarkerEdgeColor = 'none'; 

xlabel('Release time [s]');
ylabel('Delay [s]');
legend('Direction 0','Direction 1');

% Average delay direction 0
ax = gca;
Y = statistics.averageDelay_0;
line(ax.XLim, [Y Y], 'Color','black','LineStyle',':','LineWidth',2);
% Average delay direction 1
Y = statistics.averageDelay_1;
line(ax.XLim, [Y Y], 'Color', 'blue','LineStyle',':','LineWidth',2);

text = ['Delays '];
try
    if ~isempty(settings.general.caseName)
        text = [text settings.general.caseName];
    end
    text = [text settings.general.subName];
catch
    text = [text 'ERROR FOR CASE NAME'];
end
title(text);

hold off




end