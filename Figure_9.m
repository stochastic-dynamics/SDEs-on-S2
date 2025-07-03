% MATLAB code for generating Figure 9 of the manuscript
% Load the data (.mat file) for each integration scheme generated from Fig. 8
clc
clear
close all
%
%%
load Benchmark
Xbench = xx;
%
load g-Weak 1.0
Xw1 = xx;
%
load g-Weak 2.0
Xw2 = xx;
%
load g-Weak 3.0
Xw3 = xx;
%
compute_pdf_L2_error( Xbench, Xw1 )
%
compute_pdf_L2_error( Xbench, Xw2 )
%
compute_pdf_L2_error( Xbench, Xw3 )
%
% End