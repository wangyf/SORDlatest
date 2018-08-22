L = [8000,4000];
dx = 50;


nx = floor(L(1)/dx)+1;
ny= floor(L(2)/dx)+1;
nt = 801;
dt = 100/12500.;

file = 'tn';
fid=fopen(file,'r');
data = fread(fid,'single');
tn3 = reshape(data,nx,ny,nt);
fclose(fid);

it = 101;
figure(2)
clf
pcolor(tn3(:,:,1));
shading flat
colorbar

