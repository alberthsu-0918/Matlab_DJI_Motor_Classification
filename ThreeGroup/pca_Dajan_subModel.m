% After pca_Dajan_3group.m and group3_transfer.m are excuted, there are several wave classified as cluster3 (unkonw status)
% use pass wave file and unkonw wave file(pcik 10 unknown wave files which are more close to fail wave file)
% , this part is done manually, to build a sub-classified model.
clear all;

L_r = 4096*30;
N = 30; % 60 number of oberservation
FrameN = 30; %length(x1)/L;
L = L_r/FrameN;
group = 2;

OKDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\passModel\';
temp0=dir([OKDir,'*.wav']);
num_temp0=length(temp0);


%files below are hard to determine
minorfailDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\minorModel\';
temp1=dir([minorfailDir,'*.wav']);
num_temp1=length(temp1);

N = num_temp0 + num_temp1; 
x = zeros(L_r , N);

for i = 1:num_temp0
    filename=[OKDir,temp0(i).name];
    x(:, i) = wavread(filename, L_r);
end

for i = 1:num_temp1
    filename=[minorfailDir,temp1(i).name];
    x(:, num_temp0+i) = wavread(filename, L_r);
end

h = hann(L);
%h = 1;
x_ps = zeros(L/2, N);
x_h = zeros(L, N);

% average FFT spectrum
for j = 1:N
    for i = 1:FrameN
        x_h(:,j) = x((i-1)*L+1 : L*i, j);
        x_h(:,j) = x_h(:,j).*h;
        x_ps_temp = abs(fft(x_h(:,j)));
        x_ps_temp = x_ps_temp(1:L/2);
        x_ps(:,j) = x_ps(:,j) + x_ps_temp; 
    end
end


%normalize and get rid off offset of power spectrum

x_ps = 20*log10(x_ps);

%get rid of offese
scaleMean = zeros(L/2,1);
for i = 1:L/2
    tempXs = x_ps(i,:);
    scaleMean(i) = mean(tempXs);
    x_ps(i,:) = tempXs - mean(tempXs);
end

%normalized
for i = 1:L/2
    x_ps(i,:) = x_ps(i,:)./sqrt( sum(x_ps(i,:).*x_ps(i,:)) );
end



%% Find 10 peaks as characteristics 
%% useless
% psSum = zeros(L/2,1);
% for i = 1:N
%     psSum = psSum + x_ps(:,N);
% end
% 
% loc = peakfinder(psSum);
% x_ps_peak = [x_ps(loc,1) x_ps(loc,2) x_ps(loc,3) x_ps(loc,4) x_ps(loc,5) x_ps(loc,6) x_ps(loc,7) x_ps(loc,8) x_ps(loc,9) x_ps(loc,10) ];


%PCA Analysis
%[pcs, trans, evs] = princomp(x_ps');
%[pcs, trans, evs] = princomp(x_ps_peak');
%[trans,pcs,evs] = pca1(x_ps_peak);
% trans = real(trans);
% Below are SVD transpose
pcsN = 10;
[U , S, V] = svd(x_ps');
%trans = (pcs'*x_ps)';
trans = (V'*x_ps)';

hold on
for i = 1:num_temp0
    plot(trans(i,1:pcsN))
end

for i = 1:num_temp1
    plot(trans(num_temp0+i,1:pcsN), 'r')
end



%k mean cluster
cntPC = 2;
X = trans(:,1:cntPC);
[idx,c1] = kmeans(X,group)

% Find my cluster centers for pass and unknown group
% Group1: pass. Group2:unknown. 
% Group2: seems fail. Pick files 12 80 96 97 98 109 110 111 112. These files seems like fail 
c2 = zeros(group,cntPC);
c2(1,:) = sum(trans(1:num_temp0, 1:cntPC)) / (num_temp0);
c2(2,:) = sum(trans(num_temp0+1 : num_temp0 + num_temp1, 1:cntPC)) / (num_temp1);

% t = 1:51200/L:51200/2;

% %% SVD Test
% for i = 1:5
%     meanV(:,i) = mean(x(:,i))*ones(L,1);  
% end
% xz = x - meanV;
% %[U , S, V] = svds(xz,SensorN);
% [U , S, V] = svds(xz,1);
% xr = U*S*V';
% xr_fs = abs(fft(xr(:,1)));
% xr_fs = xr_fs(1:length(xr_fs)/2);



