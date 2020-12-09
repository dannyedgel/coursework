%{
    This file is used to conduct all necessary tasks for problem set 4 of
    the second quarter of Econ 712

    Date created:  08 Dec 2020
    Last modified: 08 Dec 2020
    Author: Danny Edgel
%}

%%% clean workspace
clc; clear

%%% set VF iteration variables
tol = 1e-3;
inckap = 0.005;                   % size of capital grid increments

%%% choose which files to run
q1 = 0;
q2 = 1;

%%% run files according to selection
if q1 == 1; ps4_q1; end
if q2 == 1; ps4_q2; end







