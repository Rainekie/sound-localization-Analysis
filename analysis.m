clear all;
close all;
clc;

%% vars
% Change the number to 1 to display specific plot
plotflags = [1, 0, 0, 0, 0, 0]; %Sphere, Mercator, Mollweide, Distance diff, Azi/Ele diff, RT diff
 filename = ["Sean.csv", "MJ.csv", "Parn.csv", "Jongho.csv",  "Hongjun.csv", "Youngjun.csv"];%"Jieun.csv", "Kangeun.csv",
%filename = ["Kangeun.csv"];

HRTFtarget = ["Generic", "3D_Based", "MIT_KEMAR"];
referencePoints = readmatrix('refPoints.csv');

%% read csv
for i = 1:length(filename)
    dat = readtable(filename(i));

    if i == 1
        resArray = readcsv_forExp(dat, HRTFtarget);
    else
        tmp = readcsv_forExp(dat, HRTFtarget);
        resArray = [resArray; tmp];
    end
end

%% plot 1/2/3. Sphere, Mercator, and Mollweide Projection
for HRTFs = 1:3
    data = resArray(:,:,HRTFs);
    if plotflags(1)
        sphereplot(data, referencePoints, HRTFtarget(HRTFs));
    elseif plotflags(2)
        hoge{HRTFs} = mercatorplot(data, referencePoints, HRTFtarget(HRTFs));
    elseif plotflags(3)
        mollweideplot(data, referencePoints, HRTFtarget(HRTFs));
    end
end

%% plot4/5/6. Distance/Azimuth and Elevation/Reaction time differences
statdata = diffplot(resArray, referencePoints, HRTFtarget, plotflags(4:6));

%% circular statistics
data1 = statdata(:,:,1);
data2 = statdata(:,:,2);
data3 = statdata(:,:,3);

% 各データに矢の種類の識別子を追加
arrow_type1 = ones(size(data1, 1), 1); % 矢の種類1
arrow_type2 = 2 * ones(size(data2, 1), 1); % 矢の種類2
arrow_type3 = 3 * ones(size(data3, 1), 1); % 矢の種類3

% データと矢の種類を結合
data_with_types1 = [data1, arrow_type1];
data_with_types2 = [data2, arrow_type2];
data_with_types3 = [data3, arrow_type3];

% 全てのデータを結合
combined_data = [data_with_types1; data_with_types2; data_with_types3];

% 大円距離（ラジアン）データと矢の種類識別子を抽出
alpha = combined_data(:, 1); % 大円距離（ラジアン）
idx = combined_data(:, 2); % 矢の種類識別子

% circ_wwtest関数の呼び出し
[pval, table] = circ_wwtest(alpha, idx);

% 結果の表示
disp('p-value:');
disp(pval);

disp('ANOVA Table:');
disp(table);

% 被験者数（例: 各ファイルの列数）
numSubjects = size(data1, 2);

% データの整形
alpha = []; % 大円距離（ラジアン）データ
idp = []; % 矢の種類識別子
idq = []; % ターゲットの種類識別子

% 矢の種類とターゲットの種類の識別子を設定
for i = 1:numSubjects
    alpha = [alpha; data1(:, i); data2(:, i); data3(:, i)];
    idp = [idp; ones(size(data1, 1), 1); 2 * ones(size(data2, 1), 1); 3 * ones(size(data3, 1), 1)];
    idq = [idq; i * ones(size(data1, 1), 1); i * ones(size(data2, 1), 1); i * ones(size(data3, 1), 1)];
end

% 相互作用効果を含めるかどうか
inter = 1;

% 要因の名前
fn = {'ArrowType', 'TargetType'};

% circ_hktest関数の呼び出し
[pval, table] = circ_hktest(alpha, idp, idq, inter, fn);

% 結果の表示
disp('p-values:');
disp(pval);

disp('ANOVA Table:');
disp(table);