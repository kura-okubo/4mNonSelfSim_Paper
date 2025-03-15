
import os
import numpy as np
import fnmatch
import warnings
import parse
from tqdm import tqdm

def SBENCHreader(Tstart,Tlength,pname,fname, fs_read=0, pretrigger=False):
	""" SBENCHreader: Read binary file recorded by SBENCH 32Ch.
	

	[How to run]
	datmat,tmat,Ch_info=SBENCHreader(Tstart, Tlength, pname, fname)
	datmat,tmat,Ch_info=SBENCHreader(Tstart, Tlength, pname, fname, fs_read=1e6, pretrigger=true)
			
	[Input]
	Tstart     ::Float:   Time to start reading (s)
	Tlength    ::Float:   Time length to read (s)
	pname      ::String:  Path name (e.g. '/Users/myaccount/data')
	fname      ::String:  File name (e.g. 'FB03-BD-01')
	                            
	[Option]
	fs_read    ::Float:   Sampling frequency of reading file, used for the sake of
	                      reading speed. Default is same as recording system (10e6 Hz).

	pretrigger ::Bool:    True if including pretrigger samples. If true, the output
	                      data includes from (Tstart - pretrigger time window) to (Tstart + Tlength -
	                      pretrigger time window). If false, output from (Tstart) to (Tstart +
	                      Tlength). Default is false.

	[Output]
	datmat  ::Array:      Data matrix (V) [CH1 CH2 ... CH32]
	tmat    ::Vector:     Time matrix (s)
	Ch_info ::Struct:     Channel meta data.

	[Input Data format]
	*Rename option is deplicated in python version. So rename the binary files as follows:
	%
	/Users/myaccount/data
	|-- FB03-BD-01.AE1.bin
	|-- FB03-BD-01.AE1_header.txt
	|-- FB03-BD-01.AE2.bin
	|-- FB03-BD-01.AE2_header.txt
	|-- FB03-BD-01.AE3.bin
	|-- FB03-BD-01.AE3_header.txt
	|-- FB03-BD-01.AE4.bin
	`-- FB03-BD-01.AE4_header.txt

	2021/04/11 Kurama Okubo
	This script is extended from TEACreader.m and SBENCHreader.m.
	"""

	# Constants associated with recording systems
	nboard = 4 # number of recording board
	Nch = 8 # number of channels per board


	#=====================#
	#===Local functions===#
	#=====================#
	def get_fname(pname, fname):
	    """
	        return file name of SBENCH pure binary and header file.
	        nboard: number of board in recording system.
	    """
	    
	    finames = np.empty([nboard,2], dtype=object)
	    for i in range(nboard):
	        finames[i] = [os.path.join(pname, '{}.AE{:d}.bin'.format(fname, i+1)),
	                        os.path.join(pname, '{}.AE{:d}_header.txt'.format(fname, i+1))]
	    return finames


	def parse_SBench_header(headername):
	    """
	        Return channel meta data
	    """
	    Ch_info =  [{} for k in range(Nch)]
	    chflag = 0
	    
	    # Parse channel parameters with reading line by line
	    with open(headername) as fi:
	        tline = fi.readline()
	        while tline:
	            if fnmatch.fnmatch(tline.strip(), "[[]Ch*[]]"):
	                chflag += 1
	                # reading channel info
	                Ch_id = parse.parse("[Ch{:d}]", tline.strip())[0]
	                tline = fi.readline()
	                while not tline in ['\n', '\r\n']:
	                    tmp = tline.split("=")
	                    valname = tmp[0].strip()
	                    Ch_info[Ch_id][valname] = tmp[1]
	                    
	                    if valname in ['Name','XUnit','YUnit','Description']:
	                        Ch_info[Ch_id][valname] = tmp[1].strip()
	                    else:
	                        Ch_info[Ch_id][valname] = float(tmp[1])
	                        
	                    tline = fi.readline()

	            elif 'TSSamplerate' in tline:
	                temp = tline.split("=")
	                Samplerate = float(temp[1])
	                if Samplerate != 10e6:
	                    warnings.warn("Sampling frequency if fixed at 10MHz in BIAX board, so Samplerate is overwritten from {}Hz to 10MHz.".format(Samplerate))
	                    Samplerate = 10e6

	            elif 'Pretrigger' in tline:
	                temp = tline.split("=");
	                Pretrig_sample = int(temp[1]);
	            
	            tline = fi.readline()

	    if chflag != Nch:
	        print(chflag)
	        raise ValueError(headername+' does not contain 8 channel meta data.');  

	    for i in range(Nch):
	        Ch_info[i]["Samplerate"] = Samplerate
	        Ch_info[i]["Pretrig_sample"] = Pretrig_sample

	    return Ch_info


	#==================#
	#===Main process===#
	#==================#


	# ===== Check existance of binary and header files =====
	finames = get_fname(pname, fname)

	for i in range(finames.shape[0]):
	    binfile     = finames[i, 0]
	    headerfile  = finames[i, 1]
	    
	    if not os.path.isfile(binfile):
	        raise OSError("binary file: {} not found.".format(binfile))
	    elif not os.path.isfile(headerfile):
	        raise OSError("header file: {} not found.".format(headerfile))


	# ===== Reading header files =====
	Ch_info = [] # append dictionary into this list
	all_Samplerate = []
	all_Pretrig_sample = []

	for i in range(nboard):
	    headerfile  = finames[i, 1];
	    # Read record condition from header file
	    Ch_info_tmp = parse_SBench_header(headerfile);
	    # Store Channel info into global struct
	    Ch_info = np.hstack((Ch_info, Ch_info_tmp))

	for dtmp in Ch_info:
	    all_Samplerate.append(dtmp["Samplerate"])
	    all_Pretrig_sample.append(dtmp["Pretrig_sample"])
 
	# Check if Samplerate and Pretrigger are same for all Channels
	assert all_Samplerate.count(all_Samplerate[0]) == len(all_Samplerate)
	assert all_Pretrig_sample.count(all_Pretrig_sample[0]) == len(all_Pretrig_sample)

	Samplerate = all_Samplerate[0];         # Store global samplerate
	Pretrig_sample = all_Pretrig_sample[0]; # Store global pretrigger

	if fs_read == 0:
	    fs_read = Samplerate # set default fs_read if not specified. 

	# ===== Reading data from binary =====
	if fs_read > Samplerate:
	    raise ValueError("fs_read should be smaller than the sampling rate of recording {} [Hz].".format(Samplerate))

	rfs = np.round(Samplerate/fs_read); # ratio to determine the skipping bytes.

	if Samplerate%fs_read != 0:
	    fs_read = Samplerate/rfs; # update fs_read with rounding
	    warnings.warn("Sampling frequency of reading binary is rounded as {} [Hz]".format(fs_read))

	Ndat   = np.int64(round(Tlength*fs_read))
	datmat = np.zeros((Ndat, 4*Nch))
	tmpmat = np.zeros((Ndat, Nch))
	tmat  = np.arange(Ndat)/fs_read;

	skip = (rfs-1) * 2 * Nch; #[bytes]: skip (samples  * 2bytes * channels) 
	# print(skip)

	if pretrigger:
	    # Note: offset is in bytes, so 2 bytes (16bits) * 8 channel
	    offset=np.int64(np.round(Tstart*Samplerate)*2*Nch); # offset to the beginning of pretrigger; (Tstart+Pretrigger-Pretrigger)
	    tmat = tmat - float(Pretrig_sample)/Samplerate; # shift time vector by Pretrig_sample
	else:
	    # offset to the Tstart; (Tstart+Pretrigger)
	    offset=np.int64((np.round(Tstart*Samplerate)+Pretrig_sample)*2*Nch);

	for i in tqdm(range(nboard)):

	    binfile = finames[i, 0];

	    # Read binary file
	    with open(binfile, mode='rb') as file:
	        if rfs == 1:
	            # no skipping
	            tmpmat = np.fromfile(file, dtype=np.int16, offset=offset, count=Nch*Ndat) # read binary at one time.
	        else:
	            offset_pos = offset
	            # print("skipping mode")
	            tmpmat_seg_init = np.fromfile(file, dtype=np.int16, offset=offset, count=Nch) # move to the init position
	            tmpmat_seg = [np.fromfile(file, dtype=np.int16, offset=int(skip), count=Nch) for x in range(Ndat-1)]
	            tmpmat = np.hstack((tmpmat_seg_init, np.hstack(tmpmat_seg)))

	    tmpmat_reshape = tmpmat.reshape(-1, Nch).T # reshape for data to be associated with each sensor

	    # Store data to datmat
	    for n in range(Nch):
	        globe_Chid = i*Nch + n;
	        # Compute gain from sensor V max range
	        gain = (Ch_info[globe_Chid]["MaxRange"] - Ch_info[globe_Chid]["MinRange"] ) / 2
	        datmat[:, globe_Chid] = tmpmat_reshape[n,:]*gain/2**15; # for new system with 16bit board

	    # print('Data Read Progress: {}%\n'.format(25*(i+1)))


	return (datmat,tmat,Ch_info)
