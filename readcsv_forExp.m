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

% make separate tables based on the HRTF model
for i = 1:length(HRTFtarget)
    rows = strcmp(dat.HRTFType, HRTFtarget(i));
    tmp = dat(rows, :);
    resArray(:,:,i) = [tmp.StimuliPos tmp.StimuliNum tmp.Player_sAzimuth tmp.Player_sElevation tmp.ReactionTime];
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

