clear all
nx = 101;
ny = 101;


fid=fopen('test/out/ts0','rb');
data = fread(fid,'single')/1e6;
tsc = reshape(data,nx,ny);
fclose(fid);

figure(1)
pcolor(tsc');
shading flat
colorbar;

fid=fopen('strdrop/out/ts0','rb');
data = fread(fid,'single')/1e6;
tsc = reshape(data,nx,ny);
fclose(fid);

figure(2)
pcolor(tsc');
shading flat
colorbar;