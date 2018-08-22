%clear all;
clf
%close all;


nx = 61;
ny = 61;

x = linspace(0,60,nx);
y = linspace(0,60,ny);

fid=fopen('strdrop/out/ts0','rb');
data = fread(fid,'single')/1e6;
ts0 = reshape(data,nx,ny);
fclose(fid);


fid=fopen('strdrop/out/tn0','rb');
data = fread(fid,'single')/1e6;
tn0 = reshape(data,nx,ny);
fclose(fid);

fid=fopen('strdrop/out/tse','rb');
data = fread(fid,'single')/1e6;
tse = reshape(data,nx,ny);
fclose(fid);

fid=fopen('strdrop/out/tse2','rb');
data = fread(fid,'single')/1e6;
ts2 = reshape(data,81,81);
fclose(fid);

fid=fopen('strdrop/out/tne','rb');
data = fread(fid,'single')/1e6;
tne = reshape(data,nx,ny);
fclose(fid);

fid=fopen('strdrop/in/slip.bin','rb');
data = fread(fid,'single');
su0 = reshape(data,nx,ny);
fclose(fid);

fid=fopen('strdrop/out/sue','rb');
data = fread(fid,'single');
sue = reshape(data,nx,ny);
fclose(fid);

fid=fopen('strdrop/out/sue2','rb');
data = fread(fid,'single');
su2 = reshape(data,81,81);
fclose(fid);


figure(1)
colormap(jet)
% subplot(311)
% pcolor(x,y,(su0)');
% shading flat
% colorbar
% axis equal ij
% xlim([0,60])
% ylim([0,60])
% caxis([0 2]);

subplot(311)
pcolor(x,y,log10(abs((sue-su0)./su0))');
shading flat
colorbar
axis equal ij
xlim([0,60])
ylim([0,60])
%caxis([0 2]);

stresschange = (ts0-tse)';


%% benchmark
samp=[1,1];
slip = sue'*100;
slipS = slip;
slipD = zeros(size(slipS));
rig = 3.3e10;
lam = rig;
sfac = 0.5;
[sigmaS,sigmaD,EsS,EsD] = slip2stress9(slipS,slipD,samp,rig,lam,sfac);

subplot(312)
pcolor(x,y,stresschange);
shading flat
colorbar
axis equal ij
xlim([0,60])
ylim([0,60])

subplot(313)
pcolor(x,y,((stresschange-sigmaS)));
shading flat
colorbar
axis equal ij
xlim([0,60])
ylim([0,60])

%%
figure(2)
pcolor(ts2');
shading flat
colorbar
axis equal ij
xlim([0,80])
ylim([0,80])
%caxis([0 2]);