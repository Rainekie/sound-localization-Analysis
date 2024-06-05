function out = readcsv_forExp(dat, HRTFtarget)

% add unique number at the sound postion combination
for i = 1:height(dat)
    uniqueNumArray(i) = findUniqueNum(table2array(dat(i,2)),table2array(dat(i,3)));
    uniqueStimuliNumArray(i) = findUniqueStimuliNumArray(cell2mat(table2array(dat(i,4))));
    % hoge = findUniqueStimuliNumArray(dat(i,5))
end

tmp = table(uniqueNumArray', uniqueStimuliNumArray');
tmp = renamevars(tmp,"Var1","StimuliPos");
tmp = renamevars(tmp,"Var2","StimuliNum");

dat = [dat tmp];
dat = movevars(dat, "StimuliPos", 'After', "TargetElevation");
dat = movevars(dat, "StimuliNum", 'After', "Stimuli");
clear tmp

for i = 1:height(dat)
    if strcmp(cell2mat(table2array(dat(i,7))), 'DL_Based')
        dat{i,7} = {'MIT_KEMAR'};
    end
end

% % make separate tables based on the HRTF model
% for i = 1:length(HRTFtarget)
%     rows = strcmp(dat.HRTFType, HRTFtarget(i));
%     tmp = dat(rows, :);
%     for j = 1:height(tmp)
%         disp(j)
%         switch table2array(tmp(j,1))
%             case 1
%                 disp(j)
%                 disp(tmp(j,:))
%                 tmp(j,:) = [];
%                 continue
%             case 85
%                 disp(j)
%                 disp(tmp(j,:))
%                 tmp(j,:) = [];
%                 continue
%             case 169
%                 disp(j)
%                 disp(tmp(j,:))
%                 tmp(j,:) = [];
%                 continue
%             case 253
%                 disp(j)
%                 disp(tmp(j,:))
%                 tmp(j,:) = [];
%                 continue
%         end
% 
%         if table2array(tmp(j,9)) > 180
%             tmp(j,9) = tmp(j,9) - 360;
%         end
% 
%         tmp{j,10} = tmp{j,10} * -1;
%         if table2array(tmp(j,10)) > 180
%             % disp('d')
%             tmp(j,10) = tmp(j,10) - 360;
%         end
%         tmp{j,10} = tmp{j,10} * -1;
% 
%         % if tmp{j,10} < 5 || tmp{j,10} > 10
%         %     continue
%         % end
%     end
%     resArray(:,:,i) = [tmp.StimuliPos tmp.StimuliNum tmp.Player_sAzimuth tmp.Player_sElevation tmp.ReactionTime];
% end

% Make separate tables based on the HRTF model
for i = 1:length(HRTFtarget)
    rows = strcmp(dat.HRTFType, HRTFtarget(i));
    tmp = dat(rows, :);
    rowsToDelete = []; % Initialize an array to store rows to be deleted

    for j = 1:height(tmp)
        % disp(j)
        switch table2array(tmp(j,1))
            case {1, 85, 169, 253}
                % disp(j)
                % disp(tmp(j,:))
                rowsToDelete = [rowsToDelete; j]; % Add row to the delete list
                continue
        end

        if table2array(tmp(j,9)) > 180
            tmp(j,9) = tmp(j,9) - 360;
        end

        tmp{j,10} = tmp{j,10} * -1;
        if table2array(tmp(j,10)) > 180
            tmp(j,10) = tmp(j,10) - 360;
        end
        tmp{j,10} = tmp{j,10} * -1;

        % If tmp{j,10} < 5 || tmp{j,10} > 10
        %     continue
        % end
    end

    % Delete rows after processing
    tmp(rowsToDelete, :) = [];

    % Dynamically determine the number of rows for resArray
    numRows = height(tmp);
    if i == 1
        resArray = zeros(numRows, 5, length(HRTFtarget)); % Initialize resArray with correct size
    end

    % Store results
    resArray(1:numRows, :, i) = [tmp.StimuliPos tmp.StimuliNum tmp.Player_sAzimuth tmp.Player_sElevation tmp.ReactionTime];
end

out = resArray;

    function uniqueNum = findUniqueNum(dat1, dat2)
        switch dat1
            case 0
                if ~dat2
                    uniqueNum = 1;
                else
                    uniqueNum = 2;
                end
            case 30
                if ~dat2
                    uniqueNum = 3;
                else
                    uniqueNum = 4;
                end
            case 60
                if ~dat2
                    uniqueNum = 5;
                else
                    uniqueNum = 6;
                end
            case 90
                if ~dat2
                    uniqueNum = 7;
                else
                    uniqueNum = 8;
                end
            case 120
                if ~dat2
                    uniqueNum = 9;
                else
                    uniqueNum = 10;
                end
            case 150
                if ~dat2
                    uniqueNum = 11;
                else
                    uniqueNum = 12;
                end
            case 180
                if ~dat2
                    uniqueNum = 13;
                else
                    uniqueNum = 14;
                end
            otherwise
                uniqueNum = 0;
        end
    end

    function uniqueNum1 = findUniqueStimuliNumArray(str)
        if strcmp(str, 'Footstep')
            uniqueNum1 = 0;
        elseif strcmp(str, 'PinkNoise')
            uniqueNum1 = 1;
        else
            uniqueNum1 = 9;
        end

    end

end
