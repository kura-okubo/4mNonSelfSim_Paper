function [datmat,tmat]=TEACreader(Tstart,Tlength,pname,fname)
% TEACreader: Read binary file recorded by TEAC LX-100
%
% [datmat,tmat]=TEACreader(pname,fname)
%
% [Input]
% pname: Path name
% fname: File name
% Tstart: Time to start reading (s)
% Tlength: Time length to read (s)
%
% [Output]
% datmat: Data matrix (V)
%         [CH1 CH2 ...]
% tmat: Time matrix (s)

% 2011/02/02
% 2016/12/

% ===== File & Path name setting =====
if nargin<3
    [fname,pname]=uigetfile('*.hdr','Select *.hdr file');
    fname=fname(1:end-4);
end

% ===== Read record condition =====
[Ptitle Param]=textread([pname,fname,'.hdr'],'%s %s');

for n=1:length(Ptitle)
    if(strcmp(Ptitle{n},'RATE'))
        Fs=sscanf(Param{n},'%f');
    end
    if(strcmp(Ptitle{n},'NUM_SERIES'))
        Nch=sscanf(Param{n},'%u');
    end
    if(strcmp(Ptitle{n},'FILE_TYPE'))
        if(strcmp(Param{n},'INTEGER'))
            Ftype='int16';
        else
            Ftype='int32';
        end
    end
    if(strcmp(Ptitle{n},'SLOPE'))
        N2Vcoef=textscan(Param{n},'%f','delimiter',',');
    end
end
N2Vcoef=N2Vcoef{1};

% ===== Read binary file =====
fid=fopen([pname,fname,'.dat'],'r');

% --- Offset setting ---
Ndatbit=str2num(Ftype(4:5));
Ndatbyte=Ndatbit/8;

offset=round(Tstart*Fs)*Ndatbyte*Nch;
status=fseek(fid,offset,'bof');

% --- Read file ---
datmat=(fread(fid,[Nch,round(Tlength*Fs)],Ftype))';

% --- Convert number to voltage ---
for m=1:Nch
    datmat(:,m)=datmat(:,m)*N2Vcoef(m);
end

fclose(fid);

% ===== Make time matrix =====
tmat=[0:1/Fs:(length(datmat)-1)/Fs]';
