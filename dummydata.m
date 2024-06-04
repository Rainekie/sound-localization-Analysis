% サンプルデータ生成
numSubjects = 8;
numTrials = 10;
numArrowTypes = 3;

% ランダムに生成された角度データ（ラジアン）
angles = rand(numSubjects * numTrials * numArrowTypes, 1) * 2 * pi;

% 矢の種類の識別子を生成
arrowTypes = repmat((1:numArrowTypes)', numSubjects * numTrials, 1);

% Watson-Williams 検定の実行
alpha = angles; % 角度データ
idx = arrowTypes; % 矢の種類

% Watson-Williams 検定の呼び出し
[pval, table] = circ_wwtest(alpha, idx);

% 結果の表示
disp('p-value:');
disp(pval);

disp('ANOVA Table:');
disp(table);
