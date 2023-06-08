clear all
% Specify the path to the GTFS data files
path_to_gtfs_folder = '/Users/thomasvandendorpe/Dropbox/Thesis/Code/AIP5_gitkraken/GTFS_nmbs';


% Load the data into a table
agency = readtable(fullfile(path_to_gtfs_folder, 'agency.txt'));
stops = readtable(fullfile(path_to_gtfs_folder, 'stops.txt'));
routes = readtable(fullfile(path_to_gtfs_folder, 'routes.txt'));
trips = readtable(fullfile(path_to_gtfs_folder, 'trips.txt'));
stop_times = readtable(fullfile(path_to_gtfs_folder, 'stop_times.txt'));

% Find the stop IDs for Kortrijk and Oudenaarde
kortrijkStop = num2cell(stops(strcmp(stops.stop_name,'Courtrai'),:).stop_id);
oudenaardeStop = num2cell(stops(strcmp(stops.stop_name,'Oudenaarde'),:).stop_id);

% Find the routes that pass through Kortrijk and Oudenaarde
routesToFilter = routes(ismember(routes.route_id, ...
                unique(stop_times(ismember(stop_times.stop_id, [kortrijkStop; oudenaardeStop]),:).route_id)),:);

% Filter out the trips on the routes that pass through Kortrijk and Oudenaarde
tripsToKeep = trips(~ismember(trips.route_id, routesToFilter.route_id),:);

% Create a timetable for the remaining trips
stop_times_to_keep = stop_times(ismember(stop_times.trip_id, tripsToKeep.trip_id),:);
stop_times_to_keep.arrival_time = datetime(stop_times_to_keep.arrival_time,'InputFormat','HH:mm:ss');
stop_times_to_keep.departure_time = datetime(stop_times_to_keep.departure_time,'InputFormat','HH:mm:ss');
tt = timetable(stop_times_to_keep.trip_id, stop_times_to_keep.arrival_time, ...
               stop_times_to_keep.departure_time, stop_times_to_keep.stop_id, ...
               'VariableNames',{'TripID','ArrivalTime','DepartureTime','StopID'})


