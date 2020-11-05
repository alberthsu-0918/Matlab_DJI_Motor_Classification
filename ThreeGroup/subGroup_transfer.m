% After sub model is builded in pca_Dajan_subModel.m, use that model to
% classify wave files which are determined as unkonwn(minor fail) in pca_Dajan_3group.

%clear all;
L_r = 4096*30;
N = 30; % 60 number of oberservation
FrameN = 30; %length(x1)/L;
L = L_r/FrameN;
group = 2;

OKDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\pass_asMinor1\';  % for L_r = 4096*30, FrameN = 30;
%OKDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\pass_asMinor_noavg\';  % for L_r = 4096, FrameN = 1, no average;
temp0=dir([OKDir,'*.wav']);
num_temp0=length(temp0);

%files below are hard to determine
minorfailDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\fail_asMinor1\'; % for L_r = 4096*30, FrameN = 30;
%minorfailDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\fail_asMinor_noavg\'; %for L_r = 4096*1, FrameN = 1, no average;
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

% using pre-built pca transfer model's axis center
% Therefore new data need to minus pre-built offset, scaleMean.
x_ps = 20*log10(x_ps);
offset = repmat(scaleMean, 1, N);
x_ps = x_ps - offset;

%normalized
for i = 1:L/2
    x_ps(i,:) = x_ps(i,:)./sqrt( sum(x_ps(i,:).*x_ps(i,:)) );
end


%PCA transpose
pcsN = 10;
%transT = pcs'*x_ps;
transT = (V'*x_ps);
hold on
for i = 1:num_temp0
    plot(transT(1:pcsN,i))
end
for i = num_temp0+1:num_temp0+num_temp1
    plot(transT(1:pcsN,i), 'r')
end

%clustering
PCcnt = 2;
pass = c2(1,:);
fail = c2(2,:);
cluster = zeros(N,1);
for i =1:N
    transTr = transT(1:PCcnt,i)';
    distP = transTr - pass;
    distP = sum(distP.*distP);
    distF = transTr - fail;
    distF = sum(distF.*distF);
    if(distP<distF)
        cluster(i) = 1;
    else
        cluster(i) = 2;
    end
end
cluster';

X = transT(1:2,:);
[idx,c1] = kmeans(X',2);
idx';
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



