function F=readsac(files)
% F=readsac('files')
%
% Reads a SAC-file or a set of SAC-files.
% 'files' corresponds either to strictly one SAC-file,
%    either to a set of SAC-files with the logical expression *.
% F is a structure with a size equal to the number of SAC-files read.
%    It contains ALL the words/fields of the header and the trace in
%    the field "trace".
%
% EXAMPLE
% F=readsac('2004.02.23-17.31.00.ABH.SHZ.SAC')
%   reads the file 2004.02.23-17.31.00.ABH.SHZ.SAC and saves it into
%   the structure F (dimension: 1*1).
% F=readsac('Data/*BHZ*.SAC')
%   reads all the files of the form *BHZ*.SAC in the Data directory.
%   The size of F equals to the number of these files.
%
% Program by Jérôme Vergne
% Rewritten by György Hetényi
% 10 Jan 2005


%------------------------------------------------------------------
% By default, the signals are read in mode 'r'
modlect='r';

% Determine the type of "files"
rep_files=fileparts(files);
list_files=dir(files);

if length(list_files)==0
	disp('"readsac": File(s) do not exist');
	F.tau=-1;
else

for ifile=1:length(list_files)

	nomfich=fullfile(rep_files,list_files(ifile).name);

	clear dsac;

	% Read signals
	% ------------
	% Following the signal type, reading procedure can be different:
	%	'l' - IEEE floating point with little-endian byte ordering
	%	'b' - IEEE floating point with big-endian byte ordering 
	
	[C,MAXSIZE,ENDIAN]=computer;
	if ENDIAN=='L' bool_l=0;
	%disp(ENDIAN)
	else bool_l=1;
    end
    fidev=fopen(nomfich,modlect,'l');
if fidev > 0

	h1=fread(fidev,70,'float');		% --------------
	h2=fread(fidev,40,'long');		% reading header
	h3=fread(fidev,192,'uchar');		% --------------

	% testing byte-order, h2(7) must! be 6.
	if h2(7)~=6
		bool_l=1;
		fclose(fidev);
		fidev=fopen(nomfich,modlect,'b');
		h1=fread(fidev,70,'float');
		h2=fread(fidev,40,'long');
		h3=fread(fidev,192,'uchar');
	end

	% reading trace
	tamp=fread(fidev,inf,'float');
	
	dsac.h1=h1;
	dsac.h2=h2;
	dsac.h3=h3;

	
	% PART 1: reading float-type parameters
	% ------------------------------------------
	%      %- comment are from original version, 
	%      % comments ares from SAC-homepage

	dsac.delta=h1(1);	%- sampling time interval
	dsac.depmin=h1(2);	% minimum value of dependent variable
	dsac.depmax=h1(3);	% maximum value of dependent variable
	dsac.scale=h1(4);	% multiplying scale factor for dependent variable (not currently used)
	dsac.odelta=h1(5);	% observed increment if different from nominal value
	
	dsac.b=h1(6);		%- begin time (décalage du 1er échantillon) (beginning value of independent variable)
	dsac.e=h1(7);		% ending value of independent variable
	dsac.o=h1(8);		%- event origin time (seconds relative to reference time)
	dsac.a=h1(9);		% first arrival time (seconds relative to reference time)
	dsac.internal10=h1(10);	% INTERNAL
	
	dsac.t0=h1(11);		%- user defined time picks or markers
	dsac.t1=h1(12);		%-      (seconds relative to reference time)
	dsac.t2=h1(13);		%-
	dsac.t3=h1(14);		%-
	dsac.t4=h1(15);		%-
	dsac.t5=h1(16);		%
	dsac.t6=h1(17);		%
	dsac.t7=h1(18);		%
	dsac.t8=h1(19);		%
	dsac.t9=h1(20);		%
	
	dsac.f=h1(21);		% fini or end of event time (seconds relative to reference time)
	dsac.resp0=h1(22);	% instrument response parameters (not currently used)
	dsac.resp1=h1(23);	%
	dsac.resp2=h1(24);	%
	dsac.resp3=h1(25);	%
	dsac.resp4=h1(26);	%
	dsac.resp5=h1(27);	%
	dsac.resp6=h1(28);	%
	dsac.resp7=h1(29);	%
	dsac.resp8=h1(30);	%
	dsac.resp9=h1(31);	%

	
	dsac.stla=h1(32);	%- station latitude (degrees, north positive)
	dsac.stlo=h1(33);	%- station longitude (degrees, east positive)
	dsac.stel=h1(34);	%- station elevation (meters)
	dsac.stdp=h1(35);	% station depth below surface (meters)(not currently used)
	
	dsac.evla=h1(36);	%- event latitude (degrees, north positive)
	dsac.evlo=h1(37);	%- event longitude (degrees, east positive)
	dsac.evel=h1(38);	% event elevation (meters)(not currently used)
	dsac.evdp=h1(39);	%- event depth below surface (meters)
	dsac.mag=h1(40);	% event magnitude

	
	dsac.user0=h1(41);	%- used defined variable storage area, floating!
	dsac.user1=h1(42);	%-
	dsac.user2=h1(43);	%-
	dsac.user3=h1(44);	%
	dsac.user4=h1(45);	%
	dsac.user5=h1(46);	%
	dsac.user6=h1(47);	%
	dsac.user7=h1(48);	%
	dsac.user8=h1(49);	%
	dsac.user9=h1(50);	%

	dsac.dist=h1(51);	%- station to event distance (km)
	dsac.az=h1(52);		%- event to station azimuth (degrees)
	dsac.baz=h1(53);	%- station to event azimuth (degrees)
	dsac.gcarc=h1(54);	%- station to event great circle arc length (degrees)
	dsac.internal55=h1(55);	% INTERNAL
	
	dsac.internal56=h1(56);	% INTERNAL
	dsac.depmen=h1(57);	% mean value of dependent variable
	dsac.cmpaz=h1(58);	%- component azimuth (degrees clockwise from north)
	dsac.cmpinc=h1(59);	%- component incidence angle (degrees from vertical)

	dsac.xminimum=h1(60);	% minimum value of X (spectral files only)
	dsac.xmaximum=h1(61);	% maximum value of X (spectral files only)
	dsac.yminimum=h1(62);	% minimum value of Y (spectral files only)
	dsac.ymaximum=h1(63);	% maximum value of Y (spectral files only)
	dsac.unused64=h1(64);	% UNUSED
	dsac.unused65=h1(65);	% UNUSED
	dsac.unused66=h1(66);	% UNUSED
	dsac.unused67=h1(67);	% UNUSED
	dsac.unused68=h1(68);	% UNUSED
	dsac.unused69=h1(69);	% UNUSED
	dsac.unused70=h1(70);	% UNUSED


	% PART 2: reading long-type parameters
	% ------------------------------------
	
				% GMT time corresponding to reference (0) time in file
	dsac.nzyear=h2(1);	%- year
	dsac.nzjday=h2(2);	%- julian day
	dsac.nzhour=h2(3);	%- hour
	dsac.nzmin=h2(4);	%- minute
	dsac.nzsec=h2(5);	% second
	dsac.nzmsec=h2(6);	% millisecond
	dsac.sec=h2(5)+h2(6)/1000;	%- full second (NOT an original SAC-attribute!)
	
	dsac.nvhdr=h2(7);	% header version number: 6!
	dsac.norid=h2(8);	% origin ID (CSS 3.0)
	dsac.nevid=h2(9);	% event ID (CSS 3.0)
	
	dsac.npts=h2(10);	%- number of points per data component

	dsac.internal81=h2(11);% INTERNAL
	dsac.nwfid=h2(12);	% waveform ID (CSS 3.0)
	dsac.nxsize=h2(13);	% spectral length (spectral files only)
	dsac.nysize=h2(14);	% spectral width (spectral files only)
	dsac.unused85=h2(15);	% UNUSED
	
	dsac.iftype=h2(16);	% type of file (required)(see SAC-page)
	dsac.idep=h2(17);	% type of dependent variable (see SAC-page)
	dsac.iztype=h2(18);	% reference time equivalence (see SAC-page)
	dsac.unused89=h2(19);	% UNUSED
	dsac.iinst=h2(20);	% type of recording instrument (not currently used)
	
	dsac.istreg=h2(21);	% station geogrpahic region (not currently used)
	dsac.ievreg=h2(22);	% event geographic location (not currently used)
	dsac.ievtyp=h2(23);	% type of event (see SAC-page)
	dsac.iqual=h2(24);	% quality of data (not currently used)(see SAC-page)
	dsac.isynth=h2(25);	% synthetic data flag (not currently used)
	
	dsac.imagtyp=h2(26);	% magnitude type
	dsac.imagsrc=h2(27);	% source of magnitude information
	dsac.unused98=h2(28);	% UNUSED
	dsac.unused99=h2(29);	% UNUSED
	dsac.unused100=h2(30);	% UNUSED
	
	dsac.unused101=h2(31);	% UNUSED
	dsac.unused102=h2(32);	% UNUSED
	dsac.unused103=h2(33);	% UNUSED
	dsac.unused104=h2(34);	% UNUSED
	dsac.unused105=h2(35);	% UNUSED
	
	dsac.leven=h2(36);	% TRUE if data is evenly spaced
	dsac.lspol=h2(37);	% TRUE if station components have positive polarity (left-hand rule)
	dsac.lovrok=h2(38);	% TRUE if it is okay to overwrite this file on disk
	dsac.lcalda=h2(39);	% TRUE if DIST,AZ,BAZ and GCARC are to be calculated form station and event coordinates
	dsac.unused110=h2(40);	% UNUSED

	
	% PART 3: reading uchar-type parameters
	% -------------------------------------

	imin1=min(find(h3(1:8)==0 | h3(1:8)==32));
	imin2=min(find(h3(9:24)==0 | h3(9:24)==32));
	
	dsac.kstnm=rm_blanc(h3(1:1+imin1-1))';	%- station name
	dsac.kevnm=rm_blanc(h3(9:9+imin2-1))';	%- event name
	
	dsac.khole=rm_blanc(h3(25:32))';	% hole identification if nuclear event
	dsac.ko=rm_blanc(h3(33:40))';	% event origin time identification
	dsac.ka=rm_blanc(h3(41:48))';	% first arrival time identification
	
	dsac.kt0=rm_blanc(h3(49:56))';	%- user defined time pick identifications, h1(11:20)
	dsac.kt1=rm_blanc(h3(57:64))';	%-
	dsac.kt2=rm_blanc(h3(65:72))';	%-
	dsac.kt3=rm_blanc(h3(73:80))';	%-
	dsac.kt4=rm_blanc(h3(81:88))';	%-
	dsac.kt5=rm_blanc(h3(89:96))';		%
	dsac.kt6=rm_blanc(h3(97:104))';		%
	dsac.kt7=rm_blanc(h3(105:112))';	%
	dsac.kt8=rm_blanc(h3(113:120))';	%
	dsac.kt9=rm_blanc(h3(121:128))';	%
	
	dsac.kf=rm_blanc(h3(129:136))';		% fini identification
	
	dsac.kuser0=rm_blanc(h3(137:144))';	%- user defined variable storage area
	dsac.kuser1=rm_blanc(h3(145:152))';	%-
	dsac.kuser2=rm_blanc(h3(153:160))';	%-
	
	dsac.kcmpnm=rm_blanc(h3(161:168))';	%- component name
	dsac.knetwk=rm_blanc(h3(169:176))';	% name of seismic network
	dsac.kdatrd=rm_blanc(h3(177:184))';	% date data was read onto computer
	
	dsac.kinst=rm_blanc(h3(185:192))';	%- generic name of recording instrument
	
	
	% reading trace
	% -------------
	
	dsac.trace=tamp;
	
% 	fclose(fidev);
    fclose('all');

else
	disp(['"readsac": file ' nomfich ' : reading error - file does not exist'])
	dsac.tau=-1;
end
	F(ifile)=dsac;
end
	

end
	
% -------------------------------------------------------------
	
function ch1=rm_blanc(ch)
 
% looks for whitespace in 'ch' and removes them
 if ischar(ch)
         ch1=ch(find(double(ch)~=32) & double(ch)~=0);
 else
         ch1=ch(find(ch~=32 & ch~=0));
         ch1=char(ch1);
 end 
