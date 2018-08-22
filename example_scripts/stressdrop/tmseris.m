clear all
clf

nt = 1251;

samp = [200,200];
nx = floor(20e3/samp(1))+1;
ny = floor(20e3/samp(2))+1;

x = linspace(0,20e3,nx);
y = linspace(0,20e3,ny);

fid=fopen('strdrop/out/suline','rb');
data = fread(fid,'single');
su = reshape(data,nx,ny,nt);
fclose(fid);

fid=fopen('strdrop/out/svline','rb');
data = fread(fid,'single');
sv = reshape(data,nx,ny,nt);
fclose(fid);


fid=fopen('strdrop/out/saline','rb');
data = fread(fid,'single');
sa = reshape(data,nx,ny,nt);
fclose(fid);


fid=fopen('strdrop/out/tsline','rb');
data = fread(fid,'single');
ts = reshape(data,nx,ny,nt);
fclose(fid);
%%

fid=fopen('strdrop2/out/suline','rb');
data = fread(fid,'single');
su2 = reshape(data,nx,ny,nt);
fclose(fid);

fid=fopen('strdrop2/out/svline','rb');
data = fread(fid,'single');
sv2 = reshape(data,nx,ny,nt);
fclose(fid);


fid=fopen('strdrop2/out/saline','rb');
data = fread(fid,'single');
sa2 = reshape(data,nx,ny,nt);
fclose(fid);


fid=fopen('strdrop2/out/tsline','rb');
data = fread(fid,'single');
ts2 = reshape(data,nx,ny,nt);
fclose(fid);

%%
ix = 51;
iy = 12;

figure(1)
subplot(411)
line1 = squeeze(su(ix,iy,:));
semilogy(line1,'r','linewidth',2)
hold on
line1 = squeeze(su2(ix,iy,:));
semilogy(line1,'b','linewidth',2)

subplot(412)
line2 = squeeze(sv(ix,iy,:));
semilogy(line2,'r','linewidth',2)
hold on
line2 = squeeze(sv2(ix,iy,:));
semilogy(line2,'b','linewidth',2)

subplot(413)
line3 = squeeze(sa(ix,iy,:));
semilogy(line3,'r','linewidth',2)
hold on
line3 = squeeze(sa2(ix,iy,:));
semilogy(line3,'b','linewidth',2)

subplot(414)
line4 = squeeze(ts(ix,iy,:));
semilogy(line4,'r','linewidth',2)
hold on
line4 = squeeze(ts2(ix,iy,:));
semilogy(line4,'b','linewidth',2)

