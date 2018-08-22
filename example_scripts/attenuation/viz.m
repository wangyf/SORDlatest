clear all

nx = 60;
ny = 60;
nt = 60;

dir = 'run/out/';
fid = fopen([dir,'attr11'],'r');
data = fread(fid,'single');
f1 = reshape(data,nx,ny,nt);
fclose(fid)

it = 40;
figure(1)
clf
% pcolor(squeeze(f1(:,:,it))');
% shading flat
% colorbar
ix = 30;
iy = 30;
tmp=squeeze(f1(ix,iy,:));
plot(tmp);

hold on
ix = 31;
iy = 30;
tmp=squeeze(f1(ix,iy,:));
plot(tmp);

ix = 30;
iy = 31;
tmp=squeeze(f1(ix,iy,:));
plot(tmp);

ix = 31;
iy = 31;
tmp=squeeze(f1(ix,iy,:));
plot(tmp);



% 

% nx = 61;
% ny = 61;
% nt = 60;
% fid = fopen([dir,'vx'],'r');
% data = fread(fid,'single');
% vx = reshape(data,nx,ny);
% fclose(fid)
% 
% fid = fopen([dir,'vy'],'r');
% data = fread(fid,'single');
% vy = reshape(data,nx,ny);
% fclose(fid)
% 
% f1=sqrt(vx.^2+vy.^2)';
% 
% 
% dir = 'tmp/out/';
% fid = fopen([dir,'vx'],'r');
% data = fread(fid,'single');
% vx = reshape(data,nx,ny);
% fclose(fid)
% 
% fid = fopen([dir,'vy'],'r');
% data = fread(fid,'single');
% vy = reshape(data,nx,ny);
% fclose(fid)
% 
% f2=sqrt(vx.^2+vy.^2)';
% 
% figure(2)
% clf
% colormap(jet)
% pcolor(f1-f2);
% %pcolor(f2)
% shading flat
% colorbar
% %caxis([-1e-5,1e-5])

