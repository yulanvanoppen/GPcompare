%% Setup
clearvars
close all

X = [];
Y = [];

%% First comparision: minimal difference
load('data/X1.mat')
load('data/Y1.mat')

[p1, T1] = GPcompare(X, Y, 'yLimits', [.85 1.35], ...
                           'differenceLimit', .1721, ...
                           'effectLimit', 5.8571)

%% Second comparison: observable difference
load('data/X2.mat')
load('data/Y2.mat')

[p2, T2] = GPcompare(X, Y, 'yLimits', [.85 1.35], ...
                           'differenceLimit', .1721, ...
                           'effectLimit', 5.8571)
