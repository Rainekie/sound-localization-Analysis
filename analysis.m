clear all;
close all;
clc;

%% vars
% Change the number to 1 to display specific plot
plotflags = [0, 0, 0, 1, 0, 0]; %Sphere, Mercator, Mollweide, Distance diff, Azi/Ele diff, RT diff
filename = ["Sean.csv", "Sungjoon.csv"];

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



