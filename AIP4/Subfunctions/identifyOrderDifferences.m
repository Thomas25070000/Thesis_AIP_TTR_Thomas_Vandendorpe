function nr = identifyOrderDifferences(TT1,TT2)

trains = unique(TT1.train_id);
Ntrains = length(trains);

release1 = zeros(Ntrains,1);
release2 = zeros(Ntrains,1);
cancelled1 = zeros(Ntrains,1);
cancelled2 = zeros(Ntrains,1);

for tt = 1:length(trains)
    ev1 = TT1(find(TT1.train_id == trains(tt)),:);
    ev2 = TT2(find(TT2.train_id == trains(tt)),:);
    
    if ~ev1.cancelled(1)
        release1(tt) = ev1.adjusted_arrival(1);
    else
        cancelled1(tt) = 1;
    end
    if ~ev2.cancelled(1)
        release2(tt) = ev2.adjusted_arrival(1);
    else
        cancelled2(tt) = 1;
    end
end

Q1 = zeros(Ntrains);
Q2 = zeros(Ntrains);

nr = 0;

for ii = 1:Ntrains-1
    for jj = ii+1:Ntrains
        if ~cancelled1(ii) && ~cancelled1(jj)
            Q1(ii,jj) = release1(ii) < release1(jj);
        end
        if ~cancelled2(ii) && ~cancelled2(jj)
            Q2(ii,jj) = release2(ii) < release2(jj);
        end
        
        if Q1(ii,jj) ~= Q2(ii,jj)
            nr = nr + 1;
        end
    end
end

end