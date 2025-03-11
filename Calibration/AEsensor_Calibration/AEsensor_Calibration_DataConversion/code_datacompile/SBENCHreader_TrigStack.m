function [datmat,dattrig,tmat,Ch_info]=SBENCHreader_TrigStack(sectionlength,Nevent,pname,fname, varargin)
% SBENCHreader: Read binary file recorded by SBENCH 32Ch.
%
% Read only board 02: Ch9(data) and Ch10(sync) for the frequency response test.
% Note that the preamp id is starting from 00: thus connect cable CH08 for CH9 and CH09 for CH10. 
% Return the traces of events after removal of pretrigger samples
%
% [How to run]
%
% [datmat,dattrig,tmat,Ch_info]=SBENCHreader_TrigStack(sectionlength,Nevent,pname,fname, varargin)
%
% [Input]
% sectionlength ::Int: Number of data points with the event
% Nevent     ::Int:   Number of events to read
% pname      ::String:  Path name (e.g. '/Users/myaccount/data')
% fname      ::String:  File name (e.g. 'FB03-BD-01' with 'rename' = true
%                                    or 'data_2021-03-01_13-18-22_0' with 'rename' = false)
%
% [Option]
% fs_read    ::Float:   Sampling frequency of reading file, used for the sake of
%                       reading speed. Default is same as recording system (10e6 Hz).
% pretrigger ::Bool:    True if including pretrigger samples. If true, the output
%                       data includes from (Tstart - pretrigger time window) to (Tstart + Tlength -
%                       pretrigger time window). If false, output from (Tstart) to (Tstart +
%                       Tlength). Default is false.
%
% rename     ::Bool:    file name mode. If true, the file name should be renamed in the format as
%                       described below. If false, the original names with
%                       board ids are used. Default is true.
%
% * Options are available only when the other arguments (Tstart, Tlength, ...) are given.
% * GUI pops up when calling without arguments, where options are not
%   available (i.e. all options are set in default values).
%
% [Output]
% datmat  ::Array:      Data matrix (V) [CH9 CH10]%[CH1 CH2 ... CH32]
% tmat    ::Vector:     Time matrix (s)
% Ch_info ::Struct:     Channel meta data.
%
% [Input Data format]
% If 'rename' is true, please rename the set of binary and header files output
% from SBench as follows:
%
% /Users/myaccount/data
% |-- FB03-BD-01.AE1.bin
% |-- FB03-BD-01.AE1_header.txt
% |-- FB03-BD-01.AE2.bin
% |-- FB03-BD-01.AE2_header.txt
% |-- FB03-BD-01.AE3.bin
% |-- FB03-BD-01.AE3_header.txt
% |-- FB03-BD-01.AE4.bin
% `-- FB03-BD-01.AE4_header.txt
%
% Or if 'rename' is false, the original names with board ids are used to
% read binaries such as:
%
% /Users/myaccount/data
% |-- data_2021-03-01_13-18-22_0.bin
% |-- data_2021-03-01_13-18-22_0_binheader
% |-- data_2021-03-01_13-18-22_0_sn14792.bin
% |-- data_2021-03-01_13-18-22_0_sn14792_binheader.txt
% |-- data_2021-03-01_13-18-22_0_sn14793.bin
% |-- data_2021-03-01_13-18-22_0_sn14793_binheader.txt
% |-- data_2021-03-01_13-18-22_0_sn14794.bin
% `-- data_2021-03-01_13-18-22_0_sn14794_binheader.txt
%
% 2022/06/24 Kurama Okubo
% This script is extended from SBENCHreader.m.

%% ===== Initial parameter settings =====

Nch = 8;            % channel number per file

% read only first 2 channels (data and sync), and skip the rest of data
ft  = '2*int16';    % set data type as 'Nch*int16' to use 'skip' option of fread().

% size of datmat and dattrig
datmat_all = zeros(Nevent, sectionlength); % corresponding to observed data
dattrig_all = zeros(Nevent, sectionlength);  % corresponding to sync data

% Parse options
p = inputParser;
p.addOptional('fs_read', 0);        % sampling frequency when reading binary
p.addOptional('pretrigger', false); % true if including pretrigger samples
p.addOptional('rename', true); % true if including pretrigger samples
p.parse(varargin{:});
fs_read    = p.Results.fs_read;
pretrigger = p.Results.pretrigger;
rename     = p.Results.rename;

% Call GUI for setting Tstart and Tlength
if nargin<1
    prompt = {'Enter Tstart (s)', 'Enter Tlength (s)'};
    dlgtitle = 'SBENCHreader';
    dims = [1 35];
    definput = {'0','0.005'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    Tstart = str2double(answer(1));
    Tlength = str2double(answer(2));
end

% Call GUI for setting pname and fname
if nargin<3
    [fname,pname]=uigetfile('*.AE1.bin','Select *.AE1.bin file');
    fname=fname(1:end-8);
end

%% ===== Check existance of binary and header files =====
finames = get_fname(pname, fname, rename);

% for i = 1:4
i = 2;
    binfile     = char(finames(i, 1));
    headerfile  = char(finames(i, 2));
    
    if ~isfile(binfile)
        error([binfile, ' is missing in the project directory.']);
    elseif  ~isfile(headerfile)
        error([headerfile, ' is missing in the project directory.']);
    end
% end

%% ===== Reading header files =====
Ch_info(4*Nch) = struct();
all_Samplerate = zeros(4*Nch);
all_Pretrig_sample = zeros(4*Nch);

% for i = 1:4
i = 2;
    headerfile  = char(finames(i, 2));
    % Read record condition from header file
    Ch_info_tmp = parse_SBench_header(headerfile);
    % Store Channel info into global struct
    for n=1:Nch
        globe_Chid = (i-1)*Nch + n;
        for fn = fieldnames(Ch_info_tmp(n))'
            Ch_info(globe_Chid).(fn{1}) = Ch_info_tmp(n).(fn{1});
        end
        all_Samplerate(globe_Chid) = Ch_info_tmp(n).Samplerate;
        all_Pretrig_sample(globe_Chid) = Ch_info_tmp(n).Pretrig_sample;
        
    end
% end

Samplerate = all_Samplerate(9); %all_Samplerate(1);         % Store global samplerate
Pretrig_sample = all_Pretrig_sample(9); %all_Pretrig_sample(1); % Store global pretrigger

% if fs_read == 0; fs_read = Samplerate; end % set default fs_read if not specified. 
if fs_read == 0
    fs_read = Samplerate;
elseif fs_read ~= Samplerate
    error("fs_read should be same with the Samplerate in TrigStack.");
end

%% ===== Reading data from binary =====

if fs_read > Samplerate
    error("fs_read should be smaller than the sampling rate of recording %4.2e [Hz].", Samplerate);
end

rfs = round(Samplerate/fs_read);
if rem(Samplerate, fs_read) ~= 0
    fs_read = Samplerate/rfs; % update fs_read with rounding
    warning("Sampling frequency of reading binary is rounded as %4.2e [Hz]", fs_read);
end

% skip = (rfs - 1) * 2 * 8; %[bytes]: skip (samples * bytes * channels) in fread().
skip = 2 * 6; %[bytes]: skip CH10-15 in fread().

Ndat   = int64(sectionlength*Nevent);
tmpmat = zeros(2, Ndat);


offset = 0; % we read from the beginning of data

% for i = 1:4
i = 2;
    binfile     = char(finames(i, 1));

    % Read binary file
    fid=fopen(binfile,'r');
    fseek(fid,offset,'bof');
    skip
    tmpmat=fread(fid,[2, Ndat],ft,skip);
    fclose(fid);

    % Store data to datmat
    % save data
    
    globe_Chid_data = 9; % channel 9 for the data
    % Compute gain from sensor V max range
    globe_Chid_data
    gain_data = (Ch_info(globe_Chid_data).MaxRange - Ch_info(globe_Chid_data).MinRange)/2;
    gain_data
    tr_data = tmpmat(1,:)*gain_data/2^15; % for new system
    datmat_all = reshape(tr_data, [sectionlength, Nevent])';

    % save sync
    globe_Chid_sync = 10; % channel 9 for the data
    gain_sync = (Ch_info(globe_Chid_sync).MaxRange - Ch_info(globe_Chid_sync).MinRange)/2;
    gain_sync
    tr_sync = tmpmat(2,:)*gain_sync/2^15; % for new system
    dattrig_all = reshape(tr_sync, [sectionlength, Nevent])';

% end

% Remove pretrigger
npts = sectionlength - Pretrig_sample;
tmat  = [0:double(npts-1)]'/fs_read;
datmat = datmat_all(:, Pretrig_sample+1:end);
dattrig = dattrig_all(:, Pretrig_sample+1:end);

end

%% Local functions
function [Ch_info] = parse_SBench_header(headername)
% Return channel meta data

Ch_info(8) = struct();
chflag = 0;
% Parse channel parameters with reading line by line
fid = fopen(headername);
tline = fgetl(fid);
while ischar(tline)
    if regexp(tline, regexptranslate('wildcard', '[Ch*]'))
        chflag = chflag+1;
        % reading channel info
        Ch_id = sscanf(tline,'[Ch%d]') + 1;
        tline = fgetl(fid);
        while ~isempty(tline)
            temp = split(tline, "=");
            valname = strtrim(cell2mat(temp(1)));

            if any(strcmp({'Name','XUnit','YUnit','Description'},valname))
                Ch_info(Ch_id).(valname) = strtrim(cell2mat(temp(2)));
            else
                Ch_info(Ch_id).(valname) = str2double(cell2mat(temp(2)));
            end
            
            tline = fgetl(fid);
        end
        
    elseif contains(tline, 'Samplerate') && ~contains(tline, 'TSSamplerate')
        temp = split(tline, "=");
        Samplerate = str2double(cell2mat(temp(2)));
        if Samplerate ~= 10e6
            warning("Sampling frequency if fixed at 10MHz in BIAX board, so Samplerate is overwritten from %eHz to 10MHz.", Samplerate);
            Samplerate = 10e6;
        end
        
    elseif contains(tline, 'Pretrigger')
        temp = split(tline, "=");
        Pretrig_sample = str2num(cell2mat(temp(2)));

    end
    
    tline = fgetl(fid);

end

if chflag ~= 8
    error([headername, 'does not contain 8 channel meta data.']);  
end

for i = 1:8
    Ch_info(i).Samplerate = Samplerate;
    Ch_info(i).Pretrig_sample = Pretrig_sample;
end

end

function [finames] = get_fname(pname, fname, rename)
% return file name with option 'rename'

    finames = cell(4,2);
    
    if rename
        for i = 1:4
            finames(i, :) = {fullfile(pname, sprintf('%s.AE%d.bin', fname, i)) ...
                          fullfile(pname, sprintf('%s.AE%d_header.txt', fname, i))};
        end
        
    else

        finames(1, :) = {fullfile(pname, sprintf('%s.bin', fname)) ...
                          fullfile(pname, sprintf('%s_binheader.txt', fname))};
        finames(2, :) = {fullfile(pname, sprintf('%s_sn14792.bin', fname)) ...
                          fullfile(pname, sprintf('%s_sn14792_binheader.txt', fname))};
        finames(3, :) = {fullfile(pname, sprintf('%s_sn14793.bin', fname)) ...
                          fullfile(pname, sprintf('%s_sn14793_binheader.txt', fname))};
        finames(4, :) = {fullfile(pname, sprintf('%s_sn14794.bin', fname)) ...
                          fullfile(pname, sprintf('%s_sn14794_binheader.txt', fname))};     
    end
end
