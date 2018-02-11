function [] = write_FVCOM_mhke(mhkx,mhky,At,Cp,mhkfile)

nt = numel(mhkx);

if(exist(mhkfile)); delete(mhkfile); end;
fid = fopen(mhkfile,'w');
fprintf(fid,'NTURBINES = %d\n',nt);
for i=1:nt
  fprintf(fid,'%d %f %f %f %f\n',i,mhkx(i),mhky(i),At(i),Cp(i)); 
end;
fclose(fid);
