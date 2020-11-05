% clc;
% Dajan Test.  
% After run 'pcaTest_Dajan_3group.m' to get pcs, the transpose matrix
% Run this m file to see if the pcs can project other wavefile to a better
% 3d axis.
% Further clustering is needed
% may need a sub classify to further determine wave of cluster3.

%clear all;
L_r = 65536;
FrameN = 1; %length(x1)/L;
L = L_r/FrameN;
dwtDepth = 3;
base = 'rbio6.8'; %bior6.8 rbio6.8
features = 6;


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
for i = 1:N
    x(:,i) = x(:,i).*h;
end

if dwtDepth == 3
    [A,D] = dwt(x(:,1), base);
    L1 = length(A);
    dwtLev1 = zeros(L1,N*2);
    [AA,DA] = dwt(A(:,1), base);
    L2 = length(AA);
    dwtLev2 = zeros(L2,N*2*2);
    [AAA,DAA] = dwt(AA(:,1), base);
    L3 = length(AAA);
    dwtLev3 = zeros(L3, N*2*2*2);
    x_ps_temp = zeros(L3, N*2^dwtDepth);
    x_ps = zeros(L3/2, N*2^dwtDepth);
    for i = 1:N
        [dwtLev1(:,i),dwtLev1(:,N+i)] = dwt(x(:,i), base);
    end
    
    for i = 1:N*2
        [dwtLev2(:,i), dwtLev2(:,2*N+i)] = dwt(dwtLev1(:,i), base);
    end
    
    for i = 1:N*4
        [dwtLev3(:,i), dwtLev3(:,4*N+i)] = dwt(dwtLev2(:,i), base);
    end
    
    
    for i = 1:N
        for j = 1:2^dwtDepth
            x_ps_temp(:,N*(j-1)+i) = abs(fft(dwtLev3(:,N*(j-1)+i)));
        end
    end
end

x_ps = x_ps_temp(1:length(x_ps_temp(:,1))/2,:);
x_ps = 20*log10(x_ps);

% using pre-built pca transfer model's axis center
% Therefore new data need to minus pre-built offset, scaleMean.

offset = zeros(length(x_ps(:,1)), N*2^dwtDepth);
for i = 1:2^dwtDepth
    offset(:,1+(i-1)*N:i*N) = repmat(scaleMean(:,i), 1, N);
end
x_ps = x_ps - offset;


% normalized std for each spectrum data
for i = 1:length(x_ps(:,1))
    for j = 1:2^dwtDepth
        x_ps(i,(j-1)*N+1:j*N) = x_ps(i,(j-1)*N+1:j*N)./ sqrt(sum(x_ps(i,(j-1)*N+1:j*N).*x_ps(i,(j-1)*N+1:j*N))); 
    end
end


%PCA transpose
transT = zeros(pcsN*N*2^dwtDepth, pcsN);
for i = 1:2^dwtDepth
    transT(1+(i-1)*N:i*N,1:pcsN) = (Vt(1+pcsN*(i-1):pcsN*i, :)*x_ps( :,1+(i-1)*N:i*N ))';
end

% for i = 1:2^dwtDepth
%     [U , S, V] = svd(x_ps_t( 1+(i-1)*N:i*N, : ));
%     Vt(1+pcsN*(i-1):pcsN*i, :) = V(:,pcsN)';
%     %trans = (pcs'*x_ps)';
%     %trans(1+(i-1)*N:i*N,:) = (V'*x_ps( :,1+(i-1)*N:i*N ))';
%     trans(1+(i-1)*N:i*N,1:pcsN) = (Vt(1+pcsN*(i-1):pcsN*i, :)*x_ps( :,1+(i-1)*N:i*N ))';
% end

%test: normalized projection on each pc
% for i = 1:pcsN
% transT(i,:) = transT(i,:)./S(i,i);
% end

newCords = zeros(N,pcsN*2^dwtDepth);
for i = 1:N
    for j = 1:pcsN*2^dwtDepth
        newCords(i, 1+(j-1)*pcsN:pcsN*j) = transT(i+(j-1)*N,1:pcsN);
    end
end

hold on
if pcsN == 1
    for i = 1:num_temp0
        plot(newCords(i,:))
    end
    for i = num_temp0+1:num_temp0+num_temp1
        plot(newCords(i,:), 'r')
    end
end

%clustering
PCcnt = 2^dwtDepth;
pass = c2(1,:);
fail = c2(2,:);
unkn = c2(3,:);
cluster = zeros(N,1);
for i =1:N
    transTr(1,1:PCcnt) = newCords(i,1:PCcnt);  
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



