clear all


dx = 25;
L = 22000;
nx = round(L/dx)+1;
T = 15;
dt = dx/12500;
nt = round(T/dt)+1;

folder = 'asperity/out/';
%% bidata2
fid = fopen([folder,'svm'],'rb');
data = fread(fid,'single');
svm = reshape(data,nx,nt);
fclose(fid);

fid = fopen([folder,'sum'],'rb');
sum = fread(fid,'single');
%sum = reshape(data,nx);
fclose(fid);


fid = fopen([folder,'tsm'],'rb');
data = fread(fid,'single');
tsm = reshape(data,nx,nt);
fclose(fid);

x = (0:nx-1)*dx/1e3-11;
t = (0:nt-1)*dt;

%it = 130;
ddt = 20;
for it = 1:ddt:nt
figure(1)
clf
subplot(211)
plot(x,squeeze(svm(:,it)),'r');
xlim([-11 11]);
ylim([-1 10]);

subplot(212)
plot(x,squeeze(tsm(:,it)),'b');
xlim([-11 11]);
ylim([0e6 25e6]);

%str = [folder,'figs/',num2str(it,'%03d'),'.png'];
%print('-dpng',str);
end

% ix = 276;
% figure(2)
% clf
% plot(t,squeeze(svm(ix,:)),'r')
% xlim([-5,10])