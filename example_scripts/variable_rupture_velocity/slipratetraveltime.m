clear all
%% parameter space
dx = 25;
nx = 881;
x = (0:nx-1)*dx;
x = x-1e3;
vs = 3464;
vp = 6000;
fmax = vs/10/dx;

dt = 0.002;
%dt = 0.0008;
nt = round(15/dt)+1;
%nt = 12001;
%nt = 35001;
nrec = 1;

t = (0:nt-1)*dt;
f=((0:nt-1)/(nt*dt));
%f = (0:nt-1)*0.6*vs/(dx*nx);
omega = f*2*pi;
omax = fmax*2*pi;


%nx = 2201;
k = 2*pi*(0:nx-1)/(nx*dx);
kmax = fmax*2*pi/vp;
kmax = 0.02;

folder = 'MIRA_asp1/';

%% read in sliprate file

fid = fopen([folder,'svm'],'rb');
svm = fread(fid,'single');
svm = reshape(svm,nx,nt);
fclose(fid);

fid = fopen([folder,'sum'],'rb');
sum = fread(fid,'single');
fclose(fid);

% window = sum/max(sum);
% 
% figure(100)
% clf
% %test window function
% trup = abs(x)/0.6/vs;
% spec = zeros(nx,nt);
% sk = zeros(nx,nt);
% for i = 1:nx
%     spec(i,:) = fft(svm(i,:));
%     spec(i,:) = spec(i,:)./spec(i,1);%.*exp(sqrt(-1)*2*pi.*f*trup(i));
% end
% window = hanning(nx);
% for it = 1:nt
%     sk(:,it) = fft(conj(spec(:,it).*window));
% end
% sk = conj(sk);

%%
% figure(100)
% clf
% pcolor(k,omega,log10(abs(sk))')
% shading flat
% xlim([0 0.02]);
% ylim([0 60]);
% grid on
% set(gca,'layer','top')

%%
% figure(110)
% clf
% vsh = 1e-4;
% for i = 1:nx
%     if(sum(i) > vsh)
%         svm(i,:) = svm(i,:)/sum(i); 
%     else
%         svm(i,:) = nan;
%     end
%     a = svm(i,:);
%     a(a<vsh) = nan;
%     svm(i,:) = a;
% end
% 
% 
% ii = 0;
% for i = 1:40:nx
%     ii = ii + 1;
%     plot(t-trup(i),2*svm(i,:)+ii);
%     if(ii ==1 ) hold on;end   
% end
%plot(x,window.*abs(spec(:,100)))


vsh = 0.1;
for i = 1:nx
    if(sum(i) > vsh)
        svm(i,:) = svm(i,:)/sum(i); 
    else
        svm(i,:) = nan;
    end
    a = svm(i,:);
    a(a<vsh) = nan;
    svm(i,:) = a;
end

figure(11)
clf

colormap(jet)
pcolor(x,t,(log10(svm)'));
shading flat
colorbar

hold on

%rupture velocity
v1 = 0.9 * vs;
t1 = abs(x/v1);
plot(x,t1,'linewidth',0.5,'linestyle','-.','color','w')

% P stopping phase
% v1 = vp;
% t1 = -(x/v1)+6.5;
% plot(x,t1,'linewidth',0.5,'linestyle','-.','color',rgb('black'))
% 
% v1 = vp;
% t1 = (x/v1)+6.5;
% plot(x,t1,'linewidth',0.5,'linestyle','-.','color',rgb('black'))


% S stopping phase
% v1 =  0.9*vs;
% t1 = -(x/v1) + 8;
% plot(x,t1,'linewidth',0.5,'linestyle','-.','color',rgb('black'))
% 
% v1 =  0.9*vs;
% t1 = (x/v1) + 8;
% plot(x,t1,'linewidth',0.5,'linestyle','-.','color',rgb('black'))

% extra unknown slope
% v1 =  0.45*vs;
% t1 = abs(x/v1);
% plot(x,t1,'linewidth',0.5,'linestyle','-.','color',rgb('gray'))
% 
% v1 =  0.35*vs;
% t1 = abs(x/v1);
% plot(x,t1,'linewidth',0.5,'linestyle','-.','color',rgb('gray'))

