function [base_tt, hour_tt, complete_tt, traininfo] = generateGivenTimetableComplete_case3(settings,blocksections,runningtimes,rawdata)

traintype = rawdata.text(:,1);
direction = rawdata.num(:,1);
entrytime = rawdata.num(:,2) * 60 + rawdata.num(:,3);  % in [min]
firsttimepoint = min(entrytime) * 60;           % in [s] 
maxHour = settings.disruption.duration + settings.disruption.nr_hours_after;
track = rawdata.num(:,4);
Ntracks = size(settings.tracks,1);

% For each train per hour, we generate the timetable for each and every
% possible track. TRAININFO keeps track of all characteristics.
for tr = 1:Ntracks
    lab = ['track' int2str(tr)];
    timetable.(lab) = [];
end
tt = 0; % counter
train_id_dir02 = 0;
train_id_dir03 = 0;
train_id_dir12 = 0;
train_id_dir13 = 0;

traininfo = [];
ticounter = 0;

% Assuming full speed over all block sections
speedIC = settings.trains.speed.IC / 3.6; % convert to m/s
speedR = settings.trains.speed.R / 3.6; % convert to m/s


trainidvector = zeros(size(rawdata.num,1),1);


for tr = 1:Ntracks
    lab = ['track' int2str(tr)];
    tt = 0;
    train_id_dir02 = 0;
    train_id_dir03 = 0;
    train_id_dir12 = 0;
    train_id_dir13 = 0;
    ttcounter = 0;
    
	 % Start with IC trains in direction 12
    rows = find(ismember(traintype,{'IC','THA'}) & direction == 12);
    runningIC = runningtimes.IC12.regular(tr,:);
    for rr = 1:length(rows)
        
        time = entrytime(rows(rr))*60;
		hour = floor((time - firsttimepoint)/3600);

        train_id_dir12 = train_id_dir12+1;
        
        thistrack = track(rows(rr)) == tr;
        ttcounter = rows(rr);
        
		if hour < maxHour
            train_id = hour * 1000 + train_id_dir12*100 + 12;
            trainidvector(ttcounter) = train_id;
			for bb = 1:max(blocksections.id)
                if isequal(blocksections{bb,6},{'C'})
                    continue
                end
				tt = tt+1;
				timetable.(lab)(tt).event_id = tt;
				timetable.(lab)(tt).train_id = train_id;
				timetable.(lab)(tt).train_type = traintype{rows(rr)};
				timetable.(lab)(tt).direction = 12;
				timetable.(lab)(tt).blocksection = bb;
				timetable.(lab)(tt).arrival = time;
				time = time + runningIC(bb);
				timetable.(lab)(tt).departure = time;
				timetable.(lab)(tt).running = timetable.(lab)(tt).departure - timetable.(lab)(tt).arrival;

				timetable.(lab)(tt).thistrack = thistrack;
			end
		end
    end

	 % Start with IC trains in direction 13
    rows = find(ismember(traintype,{'IC','THA'}) & direction == 13);
    runningIC = runningtimes.IC13.regular(tr,:);
    for rr = 1:length(rows)
        
        time = entrytime(rows(rr))*60;
		hour = floor((time - firsttimepoint)/3600);

        train_id_dir13 = train_id_dir13+1;
        
        thistrack = track(rows(rr)) == tr;
        ttcounter = rows(rr);
        
		if hour < maxHour
            train_id = hour * 1000 + train_id_dir13*100 + 13;
            trainidvector(ttcounter) = train_id;
			for bb = 1:max(blocksections.id)
                if isequal(blocksections{bb,6},{'B'})
                    continue
                end
				tt = tt+1;
				timetable.(lab)(tt).event_id = tt;
				timetable.(lab)(tt).train_id = train_id;
				timetable.(lab)(tt).train_type = traintype{rows(rr)};
				timetable.(lab)(tt).direction = 13;
				timetable.(lab)(tt).blocksection = bb;
				timetable.(lab)(tt).arrival = time;
				time = time + runningIC(bb);
				timetable.(lab)(tt).departure = time;
				timetable.(lab)(tt).running = timetable.(lab)(tt).departure - timetable.(lab)(tt).arrival;

				timetable.(lab)(tt).thistrack = thistrack;
			end
		end
    end


% 	% R trains in direction 1
% 	rows = find(~ismember(traintype,{'IC','THA'}) & direction == 1);
% 	runningL = runningtimes.L10.regular;
% 	for rr = 1:length(rows)
% 		time = entrytime(rows(rr))*60;
% 		hour = floor((time - firsttimepoint)/3600);
%         
%         ttcounter = rows(rr);
% 		if hour < maxHour
% 			train_id_dir1 = train_id_dir1+1;
% 			thistrack = track(rows(rr)) == tr;
%             train_id = hour * 1000 + train_id_dir1*10 + 1;
%             trainidvector(ttcounter) = train_id;
% 			for bb = 1:max(blocksections.id)
% 				tt = tt+1;
% 				timetable.(lab)(tt).event_id = tt;
% 				timetable.(lab)(tt).train_id = train_id;
% 				timetable.(lab)(tt).train_type = 'R';
% 				timetable.(lab)(tt).direction = 1;
% 				timetable.(lab)(tt).blocksection = bb;
% 				timetable.(lab)(tt).arrival = time;
% 				time = time + runningL(bb);
% 				timetable.(lab)(tt).departure = time;
% 				timetable.(lab)(tt).running = timetable.(lab)(tt).departure - timetable.(lab)(tt).arrival;
% 				timetable.(lab)(tt).thistrack = thistrack;
% 			end
% 		end
% 	end

	%% Take the other direction!
	% Start with IC trains in direction 02
    rows = find(ismember(traintype,{'IC','THA'}) & direction == 02);
    runningIC = runningtimes.IC02.regular(tr,:);
    for rr = 1:length(rows)

        time = entrytime(rows(rr))*60;
        hour = floor((time - firsttimepoint)/3600);

        train_id_dir02 = train_id_dir02+1;

        thistrack = track(rows(rr)) == tr;
        ttcounter = rows(rr);
        
        if hour < maxHour
            train_id = hour * 1000 + train_id_dir02*100+2;
            trainidvector(ttcounter) = train_id;
            for bb = max(blocksections.id):-1:1
                if isequal(blocksections{bb,6},{'C'})
                    continue
                end
                tt = tt+1;
                timetable.(lab)(tt).event_id = tt;
                timetable.(lab)(tt).train_id = train_id;
                timetable.(lab)(tt).train_type = traintype{rows(rr)};
                timetable.(lab)(tt).direction = 02;
                timetable.(lab)(tt).blocksection = bb;
                timetable.(lab)(tt).arrival = time;
                time = time + runningIC(bb);
                timetable.(lab)(tt).departure = time;
                timetable.(lab)(tt).running = timetable.(lab)(tt).departure - timetable.(lab)(tt).arrival;
                timetable.(lab)(tt).thistrack = thistrack;
            end
        end
    end


	% Start with IC trains in direction 03
    rows = find(ismember(traintype,{'IC','THA'}) & direction == 03);
    runningIC = runningtimes.IC03.regular(tr,:);
    for rr = 1:length(rows)

        time = entrytime(rows(rr))*60;
        hour = floor((time - firsttimepoint)/3600);

        train_id_dir03 = train_id_dir03+1;

        thistrack = track(rows(rr)) == tr;
        ttcounter = rows(rr);
        
        if hour < maxHour
            train_id = hour * 1000 + train_id_dir03*100+3;
            trainidvector(ttcounter) = train_id;
            for bb = max(blocksections.id):-1:1
                if isequal(blocksections{bb,6},{'B'})
                    continue
                end
                tt = tt+1;
                timetable.(lab)(tt).event_id = tt;
                timetable.(lab)(tt).train_id = train_id;
                timetable.(lab)(tt).train_type = traintype{rows(rr)};
                timetable.(lab)(tt).direction = 03;
                timetable.(lab)(tt).blocksection = bb;
                timetable.(lab)(tt).arrival = time;
                time = time + runningIC(bb);
                timetable.(lab)(tt).departure = time;
                timetable.(lab)(tt).running = timetable.(lab)(tt).departure - timetable.(lab)(tt).arrival;
                timetable.(lab)(tt).thistrack = thistrack;
            end
        end
    end

% 
%     % R trains in direction 0
% 
%     rows = find(~ismember(traintype,{'IC','THA'}) & direction == 0);
%     runningL = runningtimes.L0.regular;
%     for rr = 1:length(rows)
%         time = entrytime(rows(rr))*60;
%         hour = floor((time - firsttimepoint)/3600);
%         
%         ttcounter = rows(rr);
%         
%         if hour < maxHour
%             train_id_dir0 = train_id_dir0+1;
%             thistrack = track(rows(rr)) == tr;
%             train_id = hour * 1000 + train_id_dir0*10;
%             trainidvector(ttcounter) = train_id;
%             for bb = max(blocksections.id):-1:1
%                 tt = tt+1;
%                 timetable.(lab)(tt).event_id = tt;
%                 timetable.(lab)(tt).train_id = train_id;
%                 timetable.(lab)(tt).train_type = 'R';
%                 timetable.(lab)(tt).direction = 0;
%                 timetable.(lab)(tt).blocksection = bb;
%                 timetable.(lab)(tt).arrival = time;
%                 time = time + runningL(bb);
%                 timetable.(lab)(tt).departure = time;
%                 timetable.(lab)(tt).running = timetable.(lab)(tt).departure - timetable.(lab)(tt).arrival;
%                 timetable.(lab)(tt).thistrack = thistrack;
% 			end
% 		end
% 	end
end



% Make sure that the arrival times are all positive
for tr = 1:Ntracks
    lab = ['track' int2str(tr)];
    
    add = min(-min([timetable.(lab).arrival]),0);
    for ee = 1:size(timetable.(lab),2)
        timetable.(lab)(ee).arrival = timetable.(lab)(ee).arrival + add;
        timetable.(lab)(ee).departure = timetable.(lab)(ee).departure + add;
    end
end


%% Calculate the blocking times
t_setup = settings.TT.blocktimes.setup;
t_release_R = settings.TT.blocktimes.afterR;
t_release_IC = settings.TT.blocktimes.afterIC;

% If there is no approach, we assume the same running time as in the first
% block section!
for tr = 1:Ntracks
    lab = ['track' int2str(tr)];
    
    for ee = 1:size(timetable.(lab),2)
        % What is the approach time?
        try
            if timetable.(lab)(ee).train_id ~= timetable.(lab)(ee-1).train_id
                % This is the first action of the train.
                t_before = timetable.(lab)(ee).departure - timetable.(lab)(ee).arrival + t_setup;
            else
                % There has already been an event
                t_before = timetable.(lab)(ee-1).departure - timetable.(lab)(ee-1).arrival + t_setup;
            end
        catch
            % This is the first event of the timetable!
            t_before = timetable.(lab)(ee).departure - timetable.(lab)(ee).arrival + t_setup;
        end
        timetable.(lab)(ee).start = timetable.(lab)(ee).arrival - t_before;

        switch timetable.(lab)(ee).train_type
            case {'IC','THA'}
                timetable.(lab)(ee).finish = timetable.(lab)(ee).departure + t_release_IC;
            case 'R'
                timetable.(lab)(ee).finish = timetable.(lab)(ee).departure + t_release_R;
        end
    end

    basic_timetable.(lab) = struct2table(timetable.(lab));
end
% This generated the basic timetable.

base_tt = basic_timetable;
hour_tt = basic_timetable;
complete_tt = basic_timetable;


%% Generate the traininfo structure
trains = unique(complete_tt.track1.train_id);
for tt = 1:length(trains)
    ticounter = ticounter + 1;
    traininfo(ticounter).id = trains(tt);
    rows = find(complete_tt.track1.train_id == trains(tt));
    traininfo(ticounter).ev = complete_tt.track1.event_id(rows);
    traininfo(ticounter).type = complete_tt.track1.train_type{rows(1)};
    traininfo(ticounter).dir = mod(trains(tt),100);
    traininfo(ticounter).hour = floor(trains(tt)/1000) + 1;
    traininfo(ticounter).entry = complete_tt.track1.arrival(rows(1));
    ttnumber = find(traininfo(ticounter).id == trainidvector);
    traininfo(ticounter).allowedtracks = find(rawdata.allowedtrack(ttnumber,:) == 1);
    
    for tr = 1:Ntracks
        lab = ['track' int2str(tr)];
        present = find(complete_tt.(lab).thistrack(rows));
        if ~isempty(present)
            traininfo(ticounter).track = tr;
        end
    end
end


end