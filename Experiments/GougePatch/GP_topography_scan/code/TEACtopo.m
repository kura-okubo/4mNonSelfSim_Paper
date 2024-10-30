function [tmat, Zmat, Vout]=TEACtopo(pname,fname)
% TEACdatreadmat: Read, convert, and output file recorded by TEAC LX-100
%
% TEACdatreadmat(pname,fname)
%
% [Input]
% pname: root directory
% fname: filename
% [Output]
% tmat: time vector
% Zmat: topograpy [micro m]
% 
% 2018/08/09 modified from TEACdatreader.m
% 2020/10/12 modified from TEACdatreadout.m
% 2023/2/23  modified from TEACdatreadout.m


% ===== Channel setting =====

% Trigger OUT1
CHtrig = 2;
% Laser displacement sensor output
CHlaser=3;
Laser_gain = 20.0; %10.0; %[um/V]

% ===== File & Path name setting =====
% [fname,pname]=uigetfile('*.hdr','Select *.hdr file');
% fname=fname(1:end-4);

% ===== Read binary file =====
[datmat,tmat]=TEACreader(0,inf,pname,fname);

% ===== Convert V to displacement =====
Zmat = datmat(:, CHlaser) * Laser_gain;
Vout = datmat(:, CHtrig);


