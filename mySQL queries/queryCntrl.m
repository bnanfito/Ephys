
function [data] = queryCntrl(ageLims)

%% Make connection to database
conn = database('ephysDatabase','','');

%Set query to execute on the database
query = ['SELECT tblanimal.age, ' ...
    '	tblanimal.rearingCond, ' ...
    '	tblexperiment.experimentId, ' ...
    '	tblexperiment.unitNr, ' ...
    '	tblexperiment.experimentNr, ' ...
    '	tblexperiment.module, ' ...
    '	tblexperiment.looperNameCond1, ' ...
    '	tblexperiment.looperNameCond2, ' ...
    '	tblexperiment.looperNameCond3, ' ...
    '	tblexperiment.looperNameCond4, ' ...
    '	tblexperiment.looperNameCond5, ' ...
    '	tblrecording.recSite, ' ...
    '	tblrecording.probeId, ' ...
    '	tblrecording.penNr ' ...
    'FROM ( ( ferretphysiology.tblanimal ' ...
    'INNER JOIN ferretphysiology.tblexperiment ' ...
    'ON tblanimal.experimentid = tblexperiment.experimentid)  ' ...
    'INNER JOIN ferretphysiology.tblrecording ' ...
    'ON tblexperiment.experimentkey = tblrecording.experimentkey)  ' ...
    'WHERE tblexperiment.abort = 0 ' ...
    '	AND tblanimal.rearingCond = ''normal'' ' ...
    '	AND tblanimal.age > ' num2str(ageLims(1)) ' ' ...
    '	AND tblanimal.age < ' num2str(ageLims(2)) ' ' ...
    '	AND tblexperiment.priorMFlag = 0 ' ...
    '	AND tblexperiment.duringMFlag = 0'];

%% Execute query and fetch results
data = fetch(conn,query);

%% Close connection to database
close(conn)

%% Clear variables
clear conn query

end