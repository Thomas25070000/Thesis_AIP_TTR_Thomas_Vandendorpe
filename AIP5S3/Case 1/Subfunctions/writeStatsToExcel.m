function writeStatsToExcel(filenameOutput, sheetOutput, settings, statistics)

% Check if the file is already open
[~,~,fileOpen] = xlsfinfo(filenameOutput);
if any(strcmp(fileOpen,'Sheet1')) || any(strcmp(fileOpen,'All'))
    error('Excel file is already open. Close the file and try again.');
end

columnnames = readtable(filenameOutput, 'Sheet', sheetOutput, 'Range', 'A1:AE1', 'VariableNamingRule','preserve');
myTable = readtable(filenameOutput, 'Sheet', sheetOutput);
versionNames= size(myTable);
newrow = length(versionNames)+1; 


% newVersion = ['v' int2str(settings.general.modelVersion)];
% if settings.general.modelVersion >= 9
%     orderLimit = settings.TT.orderlimit;
%     settings.general.subName = [settings.general.subName ' max' int2str(orderLimit)];
%     newVersion = ['v' int2str(settings.general.modelVersion) ' - max. ' int2str(orderLimit)];
% elseif settings.general.modelVersion == 8
%     newVersion = ['v' int2str(settings.general.modelVersion) ' - hour after'];
% elseif settings.general.modelVersion == 7
%     newVersion = ['v' int2str(settings.general.modelVersion) ' - basic'];
% end

   

% We now know the row that has to be filled, start filling the columns
tobewritten = [];
tobewritten{1} = int2str(newrow);

% Max. delay information
col = 2;
tobewritten{col} = int2str(statistics.maxDelay);
col = 3;
tobewritten{col} = statistics.maxDelay_HHMMSS;
col = 4;
tobewritten{col} = statistics.maxDelay_HHMMSS_02;
col = 5;
tobewritten{col} = statistics.maxDelay_HHMMSS_03;
col = 6;
tobewritten{col} = statistics.maxDelay_HHMMSS_12;
col = 7;
tobewritten{col} = statistics.maxDelay_HHMMSS_13;

% Total delay information
col = 8;
tobewritten{col} = int2str(statistics.totalDelay);
col = 9;
tobewritten{col} = statistics.totalDelay_HHMMSS;
col = 10;
tobewritten{col} = statistics.totalDelay_HHMMSS_02;
col = 11;
tobewritten{col} = statistics.totalDelay_HHMMSS_03;
col = 12;
tobewritten{col} = statistics.totalDelay_HHMMSS_12;
col = 13;
tobewritten{col} = statistics.totalDelay_HHMMSS_13;

% Average delay information
col = 14;
tobewritten{col} = int2str(statistics.averageDelay);
col = 15;
tobewritten{col} = statistics.averageDelay_HHMMSS;
col = 16;
tobewritten{col} = statistics.averageDelay_HHMMSS_02;
col = 17;
tobewritten{col} = statistics.averageDelay_HHMMSS_03;
col = 18;
tobewritten{col} = statistics.averageDelay_HHMMSS_12;
col = 19;
tobewritten{col} = statistics.averageDelay_HHMMSS_13;

% Cancellation information
col = 20;
tobewritten{col} = int2str(statistics.nrCancelled);
col = 21;
tobewritten{col} = int2str(statistics.nrCancelled_dir02);
col = 22;
tobewritten{col} = int2str(statistics.nrCancelled_dir03);
col = 23;
tobewritten{col} = int2str(statistics.nrCancelled_dir12);
col = 24;
tobewritten{col} = int2str(statistics.nrCancelled_dir13);

% % CPU time
% col = find(strcmp(columnnames,'CPU time (s)'));
% tobewritten{col} = num2str(statistics.CPUtime);


xlsRange = ['A' num2str(newrow)];
tobewritten_table = cell2table(tobewritten);
writetable(tobewritten_table, filenameOutput, 'Sheet', sheetOutput, 'Range', xlsRange, 'WriteVariableNames', false);


end