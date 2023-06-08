function stops = retrieveStationInformation_case3(parameters, param_values, settings)
% Make sure that parameters are in the right order!
% I.e., stop *seg*; min dwell; penalty; min number of stops dir 0; and dir 1

counter = 0;

stops = [];

% Which segments can have stops?
segments = {'A', 'S1', 'D', 'S2', 'B'};

for cc = 1:length(segments)
    seg = segments{cc};
    % Find the label
    label = ['stop ' seg];
    
    rows = find(strcmp(parameters,label));
    
    while ~isempty(rows)
        rr = rows(1);
        % We have a stop in this segment!
        if param_values(rr) > 0 && param_values(rr+1) > 0
            counter = counter + 1;
            stops(counter).segment = seg;
            stops(counter).bswithinseg = param_values(rr);
            stops(counter).bsID = param_values(rr);
            stops(counter).dwell = param_values(rr+1);
            stops(counter).weight_dir0 = param_values(rr+2);
            stops(counter).weight_dir1 = param_values(rr+3);
            stops(counter).fractionpax_dir0 = param_values(rr+4);
            stops(counter).fractionpax_dir1 = param_values(rr+5);
            stops(counter).minLOS_dir0 = param_values(rr+6);
            stops(counter).minLOS_dir1 = param_values(rr+7);
            % Possible tracks
            possible = [];
            for tr = 1:size(settings.tracks,1)
                if param_values(rr+7+tr) == 1
                    possible = [possible tr];
                end
            end
            stops(counter).allowedtrackstostop = possible;
        end
        rows(1) = [];
    end
end


end