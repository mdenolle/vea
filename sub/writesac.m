function []=writesac(C,filename,form)
% []=writesac(C,filename,form)
%
% Writes a SAC-file into file "filename"
% 
% "C" is a SAC-structure
% "form" defines writing mode, i.e. 'le' (default) for little endian
% (PC, Intel) or 'be' for big endian (Sun, HP)
%
% See readsac.m or SAC-page for more details on fieldnames
%
% Program by Jérôme Vergne
% Rewritten by György Hetényi
% 11 Jan 2005
% Revised 16 Feb 2005
%--------------------------------------------------------------------

% Check for lilliput
if nargin <3
   [C2,MAXSIZE,ENDIAN]=computer;
   if ENDIAN=='L' form='le';
   else form='be';
   end
end

% PART 1: writing float-type parameters------------------------------
C.h1(1)=C.delta;
C.h1(2)=C.depmin;
C.h1(3)=C.depmax;
C.h1(4)=C.scale;
C.h1(5)=C.odelta;


C.h1(6)=C.b;
C.h1(7)=C.e;
C.h1(8)=C.o;
C.h1(9)=C.a;
C.h1(10)=C.internal10;

C.h1(11)=C.t0;
C.h1(12)=C.t1;
C.h1(13)=C.t2;
C.h1(14)=C.t3;
C.h1(15)=C.t4;
C.h1(16)=C.t5;
C.h1(17)=C.t6;
C.h1(18)=C.t7;
C.h1(19)=C.t8;
C.h1(20)=C.t9;

C.h1(21)=C.f;
C.h1(22)=C.resp0;
C.h1(23)=C.resp1;
C.h1(24)=C.resp2;
C.h1(25)=C.resp3;
C.h1(26)=C.resp4;
C.h1(27)=C.resp5;
C.h1(28)=C.resp6;
C.h1(29)=C.resp7;
C.h1(30)=C.resp8;
C.h1(31)=C.resp9;

C.h1(32)=C.stla;
C.h1(33)=C.stlo;
C.h1(34)=C.stel;
C.h1(35)=C.stdp;

C.h1(36)=C.evla;
C.h1(37)=C.evlo;
C.h1(38)=C.evel;
C.h1(39)=C.evdp;
C.h1(40)=C.mag;

C.h1(41)=C.user0;
C.h1(42)=C.user1;
C.h1(43)=C.user2;
C.h1(44)=C.user3;
C.h1(45)=C.user4;
C.h1(46)=C.user5;
C.h1(47)=C.user6;
C.h1(48)=C.user7;
C.h1(49)=C.user8;
C.h1(50)=C.user9;

C.h1(51)=C.dist;
C.h1(52)=C.az;
C.h1(53)=C.baz;
C.h1(54)=C.gcarc;
C.h1(55)=C.internal55;

C.h1(56)=C.internal56;
C.h1(57)=C.depmen;
C.h1(58)=C.cmpaz;
C.h1(59)=C.cmpinc;
C.h1(60)=C.xminimum;

C.h1(61)=C.xmaximum;
C.h1(62)=C.yminimum;
C.h1(63)=C.ymaximum;
C.h1(64)=C.unused64;
C.h1(65)=C.unused65;
C.h1(66)=C.unused66;
C.h1(67)=C.unused67;
C.h1(68)=C.unused68;
C.h1(69)=C.unused69;
C.h1(70)=C.unused70;

% PART 2: writing long-type parameters-------------------------------

C.h2(1)=C.nzyear;
C.h2(2)=C.nzjday;
C.h2(3)=C.nzhour;
C.h2(4)=C.nzmin;
C.h2(5)=C.nzsec;
C.h2(6)=C.nzmsec;

C.h2(7)=C.nvhdr;
C.h2(8)=C.norid;
C.h2(9)=C.nevid;
C.h2(10)=C.npts;

C.h2(11)=C.internal81;
C.h2(12)=C.nwfid;
C.h2(13)=C.nxsize;
C.h2(14)=C.nysize;
C.h2(15)=C.unused85;

C.h2(16)=C.iftype;
C.h2(17)=C.idep;
C.h2(18)=C.iztype;
C.h2(19)=C.unused89;
C.h2(20)=C.iinst;

C.h2(21)=C.istreg;
C.h2(22)=C.ievreg;
C.h2(23)=C.ievtyp;
C.h2(24)=C.iqual;
C.h2(25)=C.isynth;

C.h2(26)=C.imagtyp;
C.h2(27)=C.imagsrc;
C.h2(28)=C.unused98;
C.h2(29)=C.unused99;
C.h2(30)=C.unused100;

C.h2(31)=C.unused101;
C.h2(32)=C.unused102;
C.h2(33)=C.unused103;
C.h2(34)=C.unused104;
C.h2(35)=C.unused105;

C.h2(36)=C.leven;
C.h2(37)=C.lspol;
C.h2(38)=C.lovrok;
C.h2(39)=C.lcalda;
C.h2(40)=C.unused110;

% PART 3: writing uchar-type parameters------------------------------

C.h3(1:8)=add_blanc(C.kstnm,8);
C.h3(9:24)=add_blanc(C.kevnm,16);

C.h3(25:32)=add_blanc(C.khole,8);
C.h3(33:40)=add_blanc(C.ko,8);
C.h3(41:48)=add_blanc(C.ka,8);

C.h3(49:56)=add_blanc(C.kt0,8);
C.h3(57:64)=add_blanc(C.kt1,8);
C.h3(65:72)=add_blanc(C.kt2,8);
C.h3(73:80)=add_blanc(C.kt3,8);
C.h3(81:88)=add_blanc(C.kt4,8);
C.h3(89:96)=add_blanc(C.kt5,8);
C.h3(97:104)=add_blanc(C.kt6,8);
C.h3(105:112)=add_blanc(C.kt7,8);
C.h3(113:120)=add_blanc(C.kt8,8);
C.h3(121:128)=add_blanc(C.kt9,8);

C.h3(129:136)=add_blanc(C.kf,8);

C.h3(137:144)=add_blanc(C.kuser0,8);
C.h3(145:152)=add_blanc(C.kuser1,8);
C.h3(153:160)=add_blanc(C.kuser2,8);

C.h3(161:168)=add_blanc(C.kcmpnm,8);
C.h3(169:176)=add_blanc(C.knetwk,8);
C.h3(177:184)=add_blanc(C.kdatrd,8);
C.h3(185:192)=add_blanc(C.kinst,8);

% PART 4: writing out---------------------------------------------

if strcmp(form,'le')
	fidev=fopen(filename,'w','l');
elseif strcmp(form,'be')
	fidev=fopen(filename,'w','b');
else
	disp('"writesac": Writing format should be either ''le'' or ''be''');
end
fwrite(fidev,C.h1,'float');
fwrite(fidev,C.h2,'long');
fwrite(fidev,C.h3,'uchar');
fwrite(fidev,C.trace,'float');
fclose(fidev);


%------------------------------------------------------------------
function ch1=add_blanc(ch,lch)

% Adds whitespace at the end of "ch" if it is not "lch" long
if length(ch)<lch
	if size(ch,1)>1
		ch=ch';
	end
	dl=lch-length(ch);
	ch=[ch char(zeros(1,dl)+32)];
end
ch1=ch;	