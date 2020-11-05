% clc;
% Dajan Test.  
% After run 'pcaTest_Dajan_3group.m' to get pcs, the transpose matrix
% Run this m file to see if the pcs can project other wavefile to a better
% 3d axis.
% Further clustering is needed
% may need a sub classify to further determine wave of cluster3.

%clear all;
L_r = 4096*30;
FrameN = 30; %length(x1)/L;
L = L_r/FrameN;


OKDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\pass\';
temp0=dir([OKDir,'*.wav']);
num_temp0=length(temp0);



%files below are fail for sure
failDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\fail\';
temp1=dir([failDir,'*.wav']);
num_temp1=length(temp1);

N = num_temp0 + num_temp1; 
x = zeros(L_r , N);

for i = 1:num_temp0
    filename=[OKDir,temp0(i).name];
    x(:, i) = wavread(filename, L_r);
end

for i = 1:num_temp1
    filename=[failDir,temp1(i).name];
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

% normalized std for each spectrum data
for i = 1:L/2
    x_ps(i,:) = x_ps(i,:)./ sqrt(sum(x_ps(i,:).*x_ps(i,:))); 
end


%PCA transpose
pcsN = 10;
%transT = pcs'*x_ps;
transT = (V'*x_ps);

%test: normalized projection on each pc
% for i = 1:pcsN
% transT(i,:) = transT(i,:)./S(i,i);
% end


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
unkn = c2(3,:);
cluster = zeros(N,1);
for i =1:N
    transTr(1,1:2) = transT(1:PCcnt,i)';  
    distP = transTr - pass;
    distP = sum(distP.*distP);
    distF = transTr - fail;
    distF = sum(distF.*distF);
    distU = transTr - unkn;
    distU = sum(distU.*distU);
    if(distP<distF && distP<distF)
        cluster(i) = 1;
    end
    if(distF<distP && distF<distU)
        cluster(i) = 2;
    end
    if(distU<distP && distU<distF)
        cluster(i) = 3;
    end
end
cluster'
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



