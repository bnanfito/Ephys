%% Automate Importing Data by Generating Code Using the Database Explorer App
% This code reproduces the data obtained using the Database Explorer app by
% connecting to a database, executing a SQL query, and importing data into the
% MATLAB(R) workspace. To use this code, add the password for connecting to the
% database in the database command.

% Auto-generated by MATLAB (R2021b) and Database Toolbox Version 10.2 on 05-Mar-2024 17:36:11

%% Make connection to database
conn = database('ephysDatabase','','');

%Set query to execute on the database
query = ['SELECT tblanimal.experimentId, ' ...
    '	tblanimal.age, ' ...
    '	tblanimal.sex, ' ...
    '	tblanimal.rearingCond, ' ...
    '	tblexperiment.unitNr, ' ...
    '	tblexperiment.experimentNr, ' ...
    '	tblmanipulation.manipDescr, ' ...
    '	tblmanipulation.manipDetail ' ...
    'FROM ( ( ferretphysiology.tblanimal ' ...
    'INNER JOIN ferretphysiology.tblexperiment ' ...
    'ON tblanimal.experimentid = tblexperiment.experimentid)  ' ...
    'INNER JOIN ferretphysiology.tblmanipulation ' ...
    'ON tblexperiment.experimentkey = tblmanipulation.experimentkey)  ' ...
    'WHERE tblanimal.rearingCond = ''normal'' ' ...
    '	AND tblmanipulation.manipDescr = ''Train+Cool'''];

%% Execute query and fetch results
data = fetch(conn,query);

%% Close connection to database
close(conn)

%% Clear variables
clear conn query

animals = unique(data.experimentId);

