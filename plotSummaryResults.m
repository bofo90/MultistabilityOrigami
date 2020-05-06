function plotSummaryResults(unitCell,extrudedUnitCell,opt)

extraName = sprintf('/kh%2.3f_kta%2.3f_ke%2.3f_kf%2.3f', opt.Khinge,opt.KtargetAngle,opt.Kedge, opt.Kface);
folderResults = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/mat', opt.saveFile, extraName)

allFiles = dir(folderResults);
directories = 0;
succesfullFiles = 0;

for ct = 1:length(allFiles)
    if allFiles(ct).isdir || strcmp(allFiles(ct).name(1:8), 'metadata')
        % skip all directories and metadata file
        directories = directories+1;
        continue;
    end

    % parse the file name to get back hinge set
    hingeSet = getHingeSet(allFiles(ct).name);
    if ~isequal(hingeSet, opt.angleConstrFinal(1).val(:,1)) && strcmp(opt.readAngFile,'off')
        continue;
    end
    extrudedUnitCell.angleConstr = [hingeSet(:), -pi*0.985 * ones(length(hingeSet), 1)];
    % load results from file
    load(strcat(folderResults,'/', allFiles(ct).name));
    succesfullFiles = succesfullFiles + 1;
    fprintf('Plot of Hinges number %d/%d\n', succesfullFiles, length(allFiles)-directories);
                


    [CM, Radios, Stdev, EhingeInt, SumIntAngles, SumExtAngles, maxStrech, minStrech] =...
                                getData(extrudedUnitCell, opt, result);
    %                     EhingeInt = startEndValues(EhingeInt, result);
    %                     CM = startEndValues(CM, result);
    %                     Radios = startEndValues(Radios, result);
    %                     Stdev = startEndValues(Stdev, result);
    %                     maxStrech = startEndValues(maxStrech, result);
    %                     minStrech = startEndValues(minStrech, result);
    Energies = [ones(length(result.E),1)*(ct-directories), result.Eedge,...
        result.Eface, result.Ehinge, result.EtargetAngle, EhingeInt, result.exfl]
    PosStad = [ones(length(result.E),1,1)*(ct-directories),...
        CM(:,:),Radios, Stdev,maxStrech, minStrech, SumIntAngles, SumExtAngles];
    Hinges = [num2str(ct-directories),',',mat2str(hingeSet')];
    AllAngles = [extrudedUnitCell.theta result.deform(end).theta]';
    AllAngles = [ones(size(AllAngles,1),1)*(ct-directories) AllAngles];
end