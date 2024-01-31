% Check DataBase

%% Make connection to database
conn = database('ephysDatabase','','');
if isempty(conn.Message)
disp('db connection successful')
else
disp('db connection UNsuccessful')
end

%Set query to execute on the database
query = ['SELECT tblexperiment.experimentKey, tblexperiment.experimentId, tblexperiment.unitNr, tblexperiment.experimentNr, tblexperiment.module, tblparams.pname, tblparams.pval ' ...
    'FROM ( ferretphysiology.tblexperiment INNER JOIN ferretphysiology.tblparams ON tblexperiment.experimentkey = tblparams.experimentkey)  ' ...
    'WHERE tblexperiment.module = ''PG'''];

%% Execute query and fetch results
data = fetch(conn,query);

%% Find expts that do or do not have entries for parameter
param = 'Leye_bit'; newVal = '1';
goodCount = 0;
badCount = 0;
keys = unique(vertcat(data.experimentKey));
for k = 1:length(keys)
    ind = data.experimentKey == keys(k);
    params = data.pname(ind);
    if ismember(param,params)
        goodCount = goodCount+1;
        goodKeys(goodCount) = keys(k);
    else
        badCount = badCount+1;
        badKeys(badCount) = keys(k);
    end
end

% fixCount = 0;
% for e = 1:badCount
%     curKey = badKeys(e);
%     query = ['SELECT tblexperiment.experimentKey, tblexperiment.experimentID, tblexperiment.unitNr, tblexperiment.experimentNr, tblexperiment.module, tblparams.paramsKey, tblparams.pname, tblparams.pval ' ...
%         'FROM (ferretphysiology.tblexperiment INNER JOIN ferretphysiology.tblparams ON tblexperiment.experimentKey = tblparams.experimentKey) ' ...
%         'WHERE tblexperiment.experimentKey = ' num2str(curKey) ' AND tblparams.pname = ''Leye_bit'' '];
% 
% %     query = ['SELECT tblexperiment.experimentKey, tblexperiment.experimentID, tblexperiment.unitNr, tblexperiment.experimentNr, tblexperiment.module, tblparams.paramsKey, tblparams.pname, tblparams.pval ' ...
% %         'FROM (ferretphysiology.tblexperiment INNER JOIN ferretphysiology.tblparams ON tblexperiment.experimentKey = tblparams.experimentKey) ' ...
% %         'WHERE tblexperiment.experimentKey = ' num2str(curKey) ' AND tblparams.pname = ' param ' '];
%     if isempty(fetch(conn,query))
%         newEntry = table(curKey,{param},{newVal},'VariableNames',{'experimentKey','pname','pval'});
%         sqlwrite(conn,'tblparams',newEntry)
%     end
%     fixCount = fixCount+1;
%     fixedKeys(fixCount) = curKey;
%     disp(['exptKey' num2str(curKey) ' fixed'])
% 
%     
% end








%% Close connection to database
close(conn)

%% Clear variables
clear conn query