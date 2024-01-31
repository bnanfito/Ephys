%Check Param
clear all
close all

params = {'Leye_bit','Reye_bit'};
mod = 'PG';
anaPath = 'Z:\EphysNew\analyzer';
count = 0;
countPG = 0;
cd(anaPath)
expts = dir;
expts = expts(vertcat(expts.isdir)); % take folders
for e = 1:height(expts)
    animal = expts(e).name;
    % ignore folders that do not begin with 'fe' or 'FE'
    if ~(   (strcmp(animal(1),'f')&&strcmp(animal(2),'e')) || ((strcmp(animal(1),'F')&&strcmp(animal(2),'E')))  )
        continue
    end
    cd(fullfile(anaPath,animal))
    files = dir;
    files = files(~vertcat(files.isdir)); %ignore folders
    for f = 1:height(files)
        fileName = files(f).name;
        % ignore files that do not begin with 'fe' or 'FE'
        if ~(   (strcmp(fileName(1),'f')&&strcmp(fileName(2),'e')) || ((strcmp(fileName(1),'F')&&strcmp(fileName(2),'E')))  )
            continue
        end
        unit = fileName(8:10);
        expt = fileName(12:14);
        load(fileName,'-mat') %load analyzer file
        count = count+1;
        out{count,1} = animal;
        out{count,2} = unit;
        out{count,3} = expt;
        if isfield(Analyzer,'modID')
            out{count,4} = Analyzer.modID;
        else
            out{count,4} = 'xx';
        end
        for i = 1:length(params)
            pInd = 0;
            for p = 1:length(Analyzer.P.param)
                if strcmp(Analyzer.P.param{1,p}{1},params{i})
                    pInd = p;
                end
            end
            if pInd == 0
                out{count,4+1+(2*(i-1))} = nan;
                out{count,4+2+(2*(i-1))} = nan;
            else
                out{count,4+1+(2*(i-1))} = Analyzer.P.param{1,pInd}{3};
                out{count,4+2+(2*(i-1))} = Analyzer.P.param{1,pInd}{4};
            end
        end

    end
    
end

varNames = {'exptID','unitNr','exptNr','module',params{1},'dbBit1',params{2},'dbBit2'};
tbl = table(vertcat(out{:,1}),vertcat(out{:,2}),vertcat(out{:,3}),vertcat(out{:,4}),vertcat(out{:,5}),vertcat(out{:,6}), ...
    vertcat(out{:,7}),vertcat(out{:,8}),'VariableNames',varNames);
tbl = tbl(sum(tbl.module==mod,2)==2,:);