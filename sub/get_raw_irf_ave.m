function [G0] = get_raw_irf_ave(dir,sta,ys,ds,ye,de,type,twin,az)

% ts: shift to apply
% ys,ds,ye,de: start and end date of stacks
%twin: length time window (times 2)
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
G0 = zeros(9,floor(twin*40));
toto =  [char(dir)  '/' char(sta) '/' int2str(ys) '_' char(ds2) '_00_00_00_' ...
         int2str(ye) '_' char(de2) '_00_00_00'];
fid1=[ char(toto) '_BHZ_BHZ.' char(type)];fid2=[ char(toto) '_BHZ_BHN.' char(type)];fid3=[ char(toto) '_BHZ_BHE.' char(type)]; 
fid4=[ char(toto) '_BHN_BHZ.' char(type)];fid5=[ char(toto) '_BHN_BHN.' char(type)];fid6=[ char(toto) '_BHN_BHE.' char(type)];  
fid7=[ char(toto) '_BHE_BHZ.' char(type)];fid8=[ char(toto) '_BHE_BHN.' char(type)];fid9=[ char(toto) '_BHE_BHE.' char(type)];
if exist(fid1,'file')~=2  || exist(fid2,'file')~=2 || exist(fid3,'file')~=2  || ...
     exist(fid4,'file')~=2  || exist(fid5,'file')~=2 || exist(fid6,'file')~=2 || ...
     exist(fid7,'file')~=2 || exist(fid8,'file')~=2 || exist(fid9,'file')~=2
 disp(fid1)
 disp('return')
 return
end
zz1=readsac(fid1);zn=readsac(fid2);ze=readsac(fid3);
nz=readsac(fid4);nn=readsac(fid5);ne=readsac(fid6);
ez=readsac(fid7);en=readsac(fid8);ee=readsac(fid9);
% if max(abs(zz1.trace))==0 || max(abs(zn.trace))==0 || max(abs(ze.trace)) == 0 || ...
%          max(abs(nz.trace))==0 || max(abs(nn.trace))==0 || max(abs(ne.trace)) == 0 || ...
%           max(abs(ez.trace))==0 || max(abs(en.trace))==0 || max(abs(ee.trace)) == 0 ;
%       return;end
L=floor(length(zz1.trace)/2);
n=1;
% remember that we computed the Green's function with Z being upward.
for k=1:twin*40
    crap = [nn.trace(L-k+1)+nn.trace(L+k-1) ne.trace(L-k+1)+ne.trace(L+k-1) -(nz.trace(L-k+1)+nz.trace(L+k-1)) ; ...
            en.trace(L-k+1)+en.trace(L+k-1) ee.trace(L-k+1)+ee.trace(L+k-1) -(ez.trace(L-k+1)+ez.trace(L+k-1)); ...
            -(zn.trace(L-k+1)+zn.trace(L+k-1)) -(ze.trace(L-k+1)+ze.trace(L+k-1)) zz1.trace(L-k+1)+zz1.trace(L+k-1)]; 
    crap = rot_nez2rtz(crap,az);
    crap = crap / 2;
    G0(:,n)=reshape(crap',9,1);n= n +1;
end
                
                
                
