function [setuptimes] = createSetupTimeMatrix_case3(timetable, blocksections, settings, traininfo)


if settings.constraints.minbuffer
    minSame = settings.TT.headway.same * 60;
    minOther = settings.TT.headway.other * 60;
else
    minSame = 0;
    minOther = 0;
end

% Number of machines == number of tracks
Nmachines = size(settings.tracks,1);
Ntrains = size(traininfo,2);
constrained_sections = [1 2 3];
non_constrained_sections = [4 5 6 7 8 9 10 11];

setuptimes.disrupted = zeros(Ntrains,Ntrains,Nmachines,Nmachines);

for mm = 1:Nmachines
    tracklab = ['track' int2str(mm)];
    TT = timetable.(tracklab);
    trains = unique(TT.train_id)';
    for ii = 1:Ntrains
        rows_ii_dir = find(TT.train_id == trains(ii));
        dir_ii = TT.direction(rows_ii_dir(1));
        for jj = 1:Ntrains
            rows_jj_dir = find(TT.train_id == trains(jj));
            dir_jj = TT.direction(rows_jj_dir(1));
            if abs(dir_ii-dir_jj)==0
                % Trains are running in the same direction. Regardless of the
                % train type (they are running at the same speed), we get the
                % same headway!).
                rows_ii = find(TT.train_id == trains(ii));
                rows_ii_closed = find(TT.train_id == trains(ii)& ismember(TT.blocksection, constrained_sections));
                t_all_ii = sum(TT.running(rows_ii));
                t_closed_ii = sum(TT.running(rows_ii_closed));
                rows_jj = find(TT.train_id == trains(jj));
                rows_jj_closed = find(TT.train_id == trains(jj)& ismember(TT.blocksection, constrained_sections));
                t_all_jj = sum(TT.running(rows_jj));
                t_closed_jj = sum(TT.running(rows_jj_closed));
                time_ii_jj = TT.finish(rows_ii) - TT.start(rows_jj);
                min_HW_ii_jj = max(time_ii_jj);
                release_jj_ii = TT.arrival(rows_jj(1))-TT.arrival(rows_ii(1));
                time_jj_ii = TT.finish(rows_jj) - TT.start(rows_ii);
                min_HW_jj_ii = max(time_jj_ii);
                release_ii_jj = TT.arrival(rows_ii(1))-TT.arrival(rows_jj(1));
               
                setuptimes.disrupted(ii,jj,mm,mm) = min_HW_ii_jj+abs(release_jj_ii)-t_closed_ii;
                setuptimes.disrupted(jj,ii,mm,mm) = min_HW_jj_ii+abs(release_ii_jj)-t_closed_jj;
%                 setuptimes.disrupted(ii,jj,mm,mm) = min_HW_ii_jj+abs(release_jj_ii)-t_all_ii;
%                 setuptimes.disrupted(jj,ii,mm,mm) = min_HW_jj_ii+abs(release_ii_jj)-t_all_jj;
            end
            default_from_sander = 12+9+6+6+17;
            if abs(dir_ii-dir_jj)==1
                % Trains are running in the same direction. Regardless of the
                % train type (they are running at the same speed), we get the
                % same headway!).
                rows_ii = find(TT.train_id == trains(ii)& ismember(TT.blocksection, constrained_sections));
                rows_jj = find(TT.train_id == trains(jj)& ismember(TT.blocksection, constrained_sections));
                t_closed_ii = sum(TT.running(rows_ii));
                t_closed_jj = sum(TT.running(rows_jj));
                time_ii_jj = TT.finish(rows_ii) - TT.start(rows_jj);
                min_HW_ii_jj = max(time_ii_jj);
                release_jj_ii = TT.arrival(rows_jj(1))-TT.arrival(rows_ii(1));
                time_jj_ii = TT.finish(rows_jj) - TT.start(rows_ii);
                min_HW_jj_ii = max(time_jj_ii);
                release_ii_jj = TT.arrival(rows_ii(1))-TT.arrival(rows_jj(1));
            
                                
                setuptimes.disrupted(ii,jj,mm,mm) = min_HW_ii_jj+release_jj_ii-t_closed_ii+default_from_sander;
                setuptimes.disrupted(jj,ii,mm,mm) = min_HW_jj_ii+release_ii_jj-t_closed_jj+default_from_sander;
                                
            end
            default_from_sander = 12 + 9; %+ 9 + 6 + 1; % did this artificially, might be some other mistake in the code (should probably change from runnning time to .finish - .start)
            if abs(dir_ii-dir_jj)==10
                if dir_ii>10
                    rows_ii_open = find(TT.train_id == trains(ii)& ismember(TT.blocksection, non_constrained_sections));
                    t_open_ii = sum(TT.running(rows_ii_open));
                    rows_jj_open = find(TT.train_id == trains(jj)& ismember(TT.blocksection, non_constrained_sections));
                    t_open_jj = sum(TT.running(rows_jj_open));
                    first_jj = rows_jj_open(1);
                    t_approach_jj = TT.arrival(first_jj)-TT.start(first_jj);
                    first_ii = rows_ii(1);
                    t_approach_ii = TT.arrival(first_ii)-TT.start(first_ii);
                    setuptimes.disrupted(ii,jj,mm,mm) = t_open_ii + t_open_jj + default_from_sander+t_approach_jj;
                    setuptimes.disrupted(jj,ii,mm,mm) = default_from_sander+t_approach_ii;
                end
                if dir_ii<10
                    rows_ii_open = find(TT.train_id == trains(ii)& ismember(TT.blocksection, non_constrained_sections));
                    t_open_ii = sum(TT.running(rows_ii_open));
                    rows_jj_open = find(TT.train_id == trains(jj)& ismember(TT.blocksection, non_constrained_sections));
                    t_open_jj = sum(TT.running(rows_jj_open));
                    first_ii = rows_ii(1);
                    t_approach_ii = TT.arrival(first_ii)-TT.start(first_ii);
                    first_jj = rows_jj_open(1);
                    t_approach_jj = TT.arrival(first_jj)-TT.start(first_jj);               
                    setuptimes.disrupted(ii,jj,mm,mm) = default_from_sander+t_approach_jj;
                    setuptimes.disrupted(jj,ii,mm,mm) = t_open_jj + t_open_ii + default_from_sander+t_approach_ii;
                end
            end
            default_from_sander = 12 + 9+ 9 + 6 + 15+22;
            if abs(dir_ii-dir_jj)==9 || abs(dir_ii-dir_jj)==11
                rows_ii = find(TT.train_id == trains(ii)& ismember(TT.blocksection, constrained_sections));
                t_closed_ii = sum(TT.running(rows_ii));
                rows_jj = find(TT.train_id == trains(jj)& ismember(TT.blocksection, constrained_sections));
                t_closed_jj = sum(TT.running(rows_jj));
                first_jj = rows_jj(1);
                t_approach_jj = TT.approachtime(first_jj);
                first_ii = rows_ii(1);
                t_approach_ii = TT.approachtime(first_ii);            
                setuptimes.disrupted(ii,jj,mm,mm) = default_from_sander + t_approach_jj;
                setuptimes.disrupted(jj,ii,mm,mm) = default_from_sander + t_approach_ii;
            end
        end
    end
end

crossing_time = 0;
constrained_without_last_12_13 = [1 2];
constrained_without_last_02_03 = [2 3];
tracklab = ['track' int2str(1)];
TT = timetable.(tracklab);
for ii = 1:Ntrains
    rows_ii_dir = find(TT.train_id == trains(ii));
    dir_ii = TT.direction(rows_ii_dir(1));
    for jj = 1:Ntrains
        rows_jj_dir = find(TT.train_id == trains(jj));
        dir_jj = TT.direction(rows_jj_dir(1));
        if abs(dir_ii-dir_jj) == 10
            setuptimes.disrupted(ii,jj,1,2) = -10000;
            setuptimes.disrupted(ii,jj,2,1) = -10000;
            setuptimes.disrupted(jj,ii,1,2) = -10000;
            setuptimes.disrupted(jj,ii,2,1) = -10000;
        end
        if dir_ii+dir_jj == 25
            rows_ii = find(TT.train_id == trains(ii)& TT.blocksection == 3);
            t_closed_ii = sum(TT.running(rows_ii));
            rows_jj = find(TT.train_id == trains(jj)& TT.blocksection == 3);
            t_closed_jj = sum(TT.running(rows_jj));
%             rows_ii = find(TT.train_id == trains(ii)& ismember(TT.blocksection, constrained_without_last_12_13));
%             t_closed_ii = sum(TT.running(rows_ii));
%             rows_jj = find(TT.train_id == trains(jj)& ismember(TT.blocksection, constrained_without_last_12_13));
%             t_closed_jj = sum(TT.running(rows_jj));
            setuptimes.disrupted(ii,jj,1,2) = crossing_time - t_closed_jj+default_from_sander;
            setuptimes.disrupted(ii,jj,2,1) = crossing_time - t_closed_jj+default_from_sander;
            setuptimes.disrupted(jj,ii,1,2) = crossing_time - t_closed_ii+default_from_sander;
            setuptimes.disrupted(jj,ii,2,1) = crossing_time - t_closed_ii+default_from_sander;
        end
        if dir_ii + dir_jj == 5
            rows_ii = find(TT.train_id == trains(ii)& ismember(TT.blocksection, constrained_without_last_02_03));
            t_closed_ii = sum(TT.running(rows_ii));
            t_approach_row_ii = find(TT.train_id == trains(ii)& TT.blocksection == 3);
            t_approach_ii = TT.approachtime(t_approach_row_ii);
            rows_jj = find(TT.train_id == trains(jj)& ismember(TT.blocksection, constrained_without_last_02_03));
            t_closed_jj = sum(TT.running(rows_jj));
            t_approach_row_jj = find(TT.train_id == trains(jj)& TT.blocksection == 3);
            t_approach_jj = TT.approachtime(t_approach_row_jj);
            setuptimes.disrupted(ii,jj,1,2) = crossing_time + t_approach_jj+default_from_sander;
            setuptimes.disrupted(ii,jj,2,1) = crossing_time + t_approach_jj+default_from_sander;
            setuptimes.disrupted(jj,ii,1,2) = crossing_time + t_approach_ii+default_from_sander;
            setuptimes.disrupted(jj,ii,2,1) = crossing_time + t_approach_ii+default_from_sander;
%             setuptimes.disrupted(ii,jj,1,2) = crossing_time + t_approach_jj - t_closed_ii;
%             setuptimes.disrupted(ii,jj,2,1) = crossing_time + t_approach_jj - t_closed_ii;
%             setuptimes.disrupted(jj,ii,1,2) = crossing_time + t_approach_ii - t_closed_jj;
%             setuptimes.disrupted(jj,ii,2,1) = crossing_time + t_approach_ii - t_closed_jj;
        end
        if abs(dir_ii-dir_jj) == 9 || abs(dir_ii-dir_jj) == 11
            rows_ii_closed = find(TT.train_id == trains(ii)& ismember(TT.blocksection, constrained_sections));
            rows_ii_closed_without_last = find(TT.train_id == trains(ii)& ismember(TT.blocksection, constrained_without_last_02_03));
            t_approach_row_ii = find(TT.train_id == trains(ii)& TT.blocksection == 3);
            t_approach_ii = TT.approachtime(t_approach_row_ii);
            t_closed_ii = sum(TT.running(rows_ii_closed));
            t_closed_without_last_ii = sum(TT.running(rows_ii_closed_without_last));
            rows_jj_closed = find(TT.train_id == trains(jj)& ismember(TT.blocksection, constrained_sections));
            rows_jj_closed_without_last = find(TT.train_id == trains(jj)& ismember(TT.blocksection, constrained_without_last_12_13));
            t_approach_row_jj = find(TT.train_id == trains(jj)& TT.blocksection == 3);
            t_approach_jj = TT.approachtime(t_approach_row_jj);
            t_closed_jj = sum(TT.running(rows_jj_closed));
            t_closed_without_last_jj = sum(TT.running(rows_jj_closed_without_last));
            if dir_ii < 10
                setuptimes.disrupted(ii,jj,1,2) = crossing_time - t_approach_jj - t_closed_without_last_ii+default_from_sander; %had to exclude an extra block here
                setuptimes.disrupted(ii,jj,2,1) = crossing_time - t_approach_jj - t_closed_without_last_ii+default_from_sander;
%                 setuptimes.disrupted(ii,jj,1,2) = crossing_time - t_closed_without_last_jj - t_closed_without_last_ii; %had to exclude an extra block here
%                 setuptimes.disrupted(ii,jj,2,1) = crossing_time - t_closed_without_last_jj - t_closed_without_last_ii;
                setuptimes.disrupted(jj,ii,1,2) = crossing_time + t_approach_jj+default_from_sander;% might have to exclude the default here
                setuptimes.disrupted(jj,ii,2,1) = crossing_time + t_approach_jj+default_from_sander;
            end
            if dir_ii > 10
                setuptimes.disrupted(ii,jj,1,2) = crossing_time + t_approach_ii+default_from_sander;
                setuptimes.disrupted(ii,jj,2,1) = crossing_time + t_approach_ii+default_from_sander;
                setuptimes.disrupted(jj,ii,1,2) = crossing_time - t_approach_ii - t_closed_without_last_jj+default_from_sander;
                setuptimes.disrupted(jj,ii,2,1) = crossing_time - t_approach_ii - t_closed_without_last_jj+default_from_sander; 
%                 setuptimes.disrupted(jj,ii,1,2) = crossing_time - t_closed_without_last_ii - t_closed_without_last_jj;
%                 setuptimes.disrupted(jj,ii,2,1) = crossing_time - t_closed_without_last_ii - t_closed_without_last_jj; 
            end
        end
    end

end



end