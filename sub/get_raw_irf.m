function G = get_raw_irf(dir,sta,ys,ds,ye,de,typ,twin,az)

% ts: shift to apply
% ys,ds,ye,de: start and end date of stacks
%twin
G=0;
ds2 = int2str(ds);
if ds<100
    ds2=['0' int2str(ds)];
    if ds < 10
        ds2=['0' ds2];
    end
end
de2=int2str(de);
if de<100
    de2=['0' int2str(de)];
    if de < 10
        de2=['0' de2];
    end
end
toto =  [char(dir)  '/' char(sta) '/' int2str(ys) '_' char(ds2) '_00_00_00_' ...
         int2str(ye) '_' char(de2) '_00_00_00' ];
fid1=[ char(toto) '_BHZ_BHZ.' char(typ)];fid2=[ char(toto) '_BHZ_BHN.' char(typ)];fid3=[ char(toto) '_BHZ_BHE.' char(typ)]; 
fid4=[ char(toto) '_BHN_BHZ.' char(typ)];fid5=[ char(toto) '_BHN_BHN.' char(typ)];fid6=[ char(toto) '_BHN_BHE.' char(typ)];  
fid7=[ char(toto) '_BHE_BHZ.' char(typ)];fid8=[ char(toto) '_BHE_BHN.' char(typ)];fid9=[ char(toto) '_BHE_BHE.' char(typ)]; 
if exist(fid1,'file')~=2  || exist(fid2,'file')~=2 || exist(fid3,'file')~=2  || ...
     exist(fid4,'file')~=2  || exist(fid5,'file')~=2 || exist(fid6,'file')~=2 || ...
     exist(fid7,'file')~=2  || exist(fid8,'file')~=2 || exist(fid9,'file')~=2 
 disp('return')
    return
end
zz1=readsac(fid1);zn=readsac(fid2);ze=readsac(fid3);
nz=readsac(fid4);nn=readsac(fid5);ne=readsac(fid6);
ez=readsac(fid7);en=readsac(fid8);ee=readsac(fid9);
L=floor(length(zz1.trace)/2);
n=1;
for k=L-twin*40+1:L+twin*40+1
    
    crap = [nn.trace(k) ne.trace(k) -nz.trace(k) ; ...
            en.trace(k) ee.trace(k) -ez.trace(k); ...
            -zn.trace(k) -ze.trace(k) zz1.trace(k)]; 
    crap = rot_nez2rtz(crap,az);
    G(1:9,n)=reshape(crap,9,1);n= n +1;
end
                
                
                