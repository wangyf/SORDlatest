clear all;
clf
%close all;


nx = 61;
ny = 61;
T = 20;
%dt = 1000/12500;
dt = 0.05;
nt = floor(T/dt)+1;

x = linspace(0,60,nx);
y = linspace(0,60,ny);

fid=fopen('strdrop/out/tsline','rb');
data = fread(fid,'single')/1e6;
ts = reshape(data,nx,ny,nt);
fclose(fid);

fid=fopen('strdrop/out/suline','rb');
data = fread(fid,'single');
sl = reshape(data,nx,ny,nt);
fclose(fid);

for it = 1:10:nt
it
figure(1)
subplot(121)
colormap(jet)
pcolor(x,y,(sl(:,:,it))');
shading flat
colorbar
axis equal ij
xlim([0,60])
ylim([0,60])
caxis([0 2]);

subplot(122)
colormap(jet)
pcolor(x,y,(ts(:,:,it)-60)');
shading flat
colorbar
axis equal ij
xlim([0,60])
ylim([0,60])
caxis([-8,8]);

pause
end