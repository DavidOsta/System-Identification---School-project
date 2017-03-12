addpath(genpath('tests'));
addpath(genpath('src'));

clear;clc; close all

% https://www.mathworks.com/help/matlab/script-based-unit-tests.html
% https://www.mathworks.com/help/matlab/ref/testsuite.html
% https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-using-classes.html

testCase1 = FrictionModelBasicTest;
res = run(testCase1);



