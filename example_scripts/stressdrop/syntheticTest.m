
clear all;
close all;
%clf

samp = [200,200];
nx = floor(20e3/samp(1))+1;
ny = floor(20e3/samp(2))+1;

x = linspace(0,20e3,nx);
y = linspace(0,20e3,ny);

slip = zeros(nx,ny);

hypo = [10e3,10e3];

d = 1; %1m
a = 2500;

for j = 1:ny
    for i = 1:nx
        ix = (i-1)*samp(1);
        iy = (j-1)*samp(2);
        r = (ix-hypo(1))^2+(iy-hypo(2))^2;
        slip(i,j) = d*exp(-r/a^2);
    end
end
   

figure(1)
%colormap(jet)
pcolor(x,y,log10(slip+1e-60)');
shading flat
colorbar;
axis ij equal
%% write into binary
fid = fopen('bin/slip.bin','wb');
fwrite(fid,slip,'single');
fclose(fid);

% fid = fopen('bin/slip.bin','rb');
% data=fread(fid,'single');
% data2 = reshape(data,nx,ny);
% fclose(fid);
%% benchmark
slipS = slip*100;
slipD = zeros(size(slipS));
rig = 3.3e10;
lam = rig;
sfac = 0.5;
[sigmaS,sigmaD,EsS,EsD] = slip2stress9(slipS,slipD,samp/1e3,rig,lam,sfac);

figure(2)
colormap(jet)
pcolor(sigmaS);
shading flat
colorbar
axis equal ij

fprintf('1- %6.4e\n',sum(sum(sigmaS.*slip))/sum(sum(slip)));
%% SORD results

fid=fopen('strdrop/out/sue','rb');
data = fread(fid,'single');
su = reshape(data,nx,ny);
fclose(fid);

fid=fopen('strdrop/out/tse','rb');
data = fread(fid,'single')/1e6;
tsc = 60-reshape(data,nx,ny);
fclose(fid);

figure(3)
pcolor(su-slip)
shading flat
colorbar
axis equal ij

figure(4)
%colormap(jet)
pcolor(abs(tsc'-sigmaS)/max(sigmaS(:)));
%pcolor(abs(tsc'-sigmaS));
shading flat
colorbar
axis equal ij

fprintf('2- %6.4e\n',sum(sum(tsc'.*su))/sum(sum(su)));
figure(5)
%colormap(jet)
pcolor(tsc');
shading flat
colorbar
axis equal ij
