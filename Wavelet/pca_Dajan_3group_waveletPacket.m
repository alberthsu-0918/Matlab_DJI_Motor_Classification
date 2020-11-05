% clc;
% Dajan Test 
% Apply PCA to data set of the spectrum with normalization of each freq bin, and get new "projection".
% Select 30 wave files, include 10 pass, 10 fail and 10 unkonw to build the
% pca model.
% The result can be used as basics of classification.

% After model is build by the m-file. Run "group3_transfer" to classify
% wave file. 

% Refer to
% http://research.med.helsinki.fi/corefacilities/proteinchem/pca_introduction_basics.pdf
% and
% https://learnche.org/pid/latent-variable-modelling/principal-component-analysis/pca-example-analysis-of-spectral-data
% Albert Hsu  2017/11/01
% ¼Ú­C ©Ò¦VµL¼Ä

clear all;

L_r = 65536;
N = 30; % 60 number of oberservation
group = 3;
FrameN = 1; %length(x1)/L;
L = L_r/FrameN;
dwtDepth = 4;
base = 'rbio5.5'; %bior6.8 rbio6.8
features = 6;

OKDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\passModel\';
temp0=dir([OKDir,'*.wav']);
num_temp0=length(temp0);



%files below are fail for sure
failDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\failModel\';
temp1=dir([failDir,'*.wav']);
num_temp1=length(temp1);



%files below are hard to determine
minorfailDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\minorModel\';
temp2=dir([minorfailDir,'*.wav']);
num_temp2=length(temp2);


N = num_temp0 + num_temp1 + num_temp2; 
x = zeros(L_r , N);

%%----------------discrete packle wavele transform-----------------
for i = 1:num_temp0
    filename=[OKDir,temp0(i).name];
    x(:, i) = wavread(filename, L_r);
end

for i = 1:num_temp1
    filename=[failDir,temp1(i).name];
    x(:, num_temp0+i) = wavread(filename, L_r);
end

for i = 1:num_temp2
    filename=[minorfailDir,temp2(i).name];
    x(:, num_temp0+ num_temp1 +i) = wavread(filename, L_r);
end

%h = hann(L/(2^dwtDepth));
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

if dwtDepth == 4
    [A,D] = dwt(x(:,1), base);
    L1 = length(A);
    dwtLev1 = zeros(L1,N*2);
    [AA,DA] = dwt(A(:,1), base);
    L2 = length(AA);
    dwtLev2 = zeros(L2,N*2*2);
    [AAA,DAA] = dwt(AA(:,1), base);
    L3 = length(AAA);
    dwtLev3 = zeros(L3, N*2*2*2);
    [AAAA, DAAA] = dwt(AAA(:,1), base);
    L4 = length(AAAA);
    dwtLev4 = zeros(L4, N*2*2*2*2);
    x_ps_temp = zeros(L4, N*2^dwtDepth);
    x_ps = zeros(L4/2, N*2^dwtDepth);
    for i = 1:N
        [dwtLev1(:,i),dwtLev1(:,N+i)] = dwt(x(:,i), base);
    end
    
    for i = 1:N*2
        [dwtLev2(:,i), dwtLev2(:,2*N+i)] = dwt(dwtLev1(:,i), base);
    end
    
    for i = 1:N*4
        [dwtLev3(:,i), dwtLev3(:,4*N+i)] = dwt(dwtLev2(:,i), base);
    end
    
    for i = 1:N*8
        [dwtLev4(:,i), dwtLev4(:,8*N+i)] = dwt(dwtLev3(:,i), base);
    end
    
    for i = 1:N
        for j = 1:2^dwtDepth
            x_ps_temp(:,N*(j-1)+i) = abs(fft(dwtLev4(:,N*(j-1)+i)));
        end
    end
end

if dwtDepth == 5
    [A,D] = dwt(x(:,1), base);
    L1 = length(A);
    dwtLev1 = zeros(L1,N*2);
    [AA,DA] = dwt(A(:,1), base);
    L2 = length(AA);
    dwtLev2 = zeros(L2,N*2*2);
    [AAA,DAA] = dwt(AA(:,1), base);
    L3 = length(AAA);
    dwtLev3 = zeros(L3, N*2*2*2);
    [AAAA, DAAA] = dwt(AAA(:,1), base);
    L4 = length(AAAA);
    dwtLev4 = zeros(L4, N*2*2*2*2);
    [AAAAA, DAAAA] = dwt(AAAA(:,1), base);
    L5 = length(AAAAA);
    dwtLev5 = zeros(L5, N*2*2*2*2*2);
    for i = 1:N
        [dwtLev1(:,i),dwtLev1(:,N+i)] = dwt(x(:,i), base);
    end
    
    for i = 1:N*2
        [dwtLev2(:,i), dwtLev2(:,2*N+i)] = dwt(dwtLev1(:,i), base);
    end
    
    for i = 1:N*4
        [dwtLev3(:,i), dwtLev3(:,4*N+i)] = dwt(dwtLev2(:,i), base);
    end
    
    for i = 1:N*8
        [dwtLev4(:,i), dwtLev4(:,8*N+i)] = dwt(dwtLev3(:,i), base);
    end
    
    for i = 1:N*16
        [dwtLev5(:,i), dwtLev5(:,16*N+i)] = dwt(dwtLev4(:,i), base);
    end
    
    PCA_feat = zeros(features*2^dwtDepth, N);
    for i = 1:N
        for j = 1:2^dwtDepth
            PCA_feat(features*(j-1)+1,i) = var(dwtLev5(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+2,i) = skewness(dwtLev5(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+3,i) = kurtosis(dwtLev5(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+4,i) = sqrt(sum((dwtLev5(:,N*(j-1)+i)).*dwtLev5(:,N*(j-1)+i))/length(dwtLev5(:,1)));
            PCA_feat(features*(j-1)+5,i) = mean(dwtLev5(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+6,i) = (max(dwtLev5(:,N*(j-1)+i))-min(dwtLev5(:,N*(j-1)+i))) / sqrt(mean([dwtLev5(:,N*(j-1)+i)].^2)); %peak 2 rms. aka crest factor
        end
    end
end

if dwtDepth == 6
    [A,D] = dwt(x(:,1), base);
    L1 = length(A);
    dwtLev1 = zeros(L1,N*2);
    [AA,DA] = dwt(A(:,1), base);
    L2 = length(AA);
    dwtLev2 = zeros(L2,N*2*2);
    [AAA,DAA] = dwt(AA(:,1), base);
    L3 = length(AAA);
    dwtLev3 = zeros(L3, N*2*2*2);
    [AAAA, DAAA] = dwt(AAA(:,1), base);
    L4 = length(AAAA);
    dwtLev4 = zeros(L4, N*2*2*2*2);
    [AAAAA, DAAAA] = dwt(AAAA(:,1), base);
    L5 = length(AAAAA);
    dwtLev5 = zeros(L5, N*2*2*2*2*2);
    [AAAAAA, DAAAAA] = dwt(AAAAA(:,1), base);
    L6 = length(AAAAAA);
    dwtLev6 = zeros(L6, N*2*2*2*2*2*2);
    for i = 1:N
        [dwtLev1(:,i),dwtLev1(:,N+i)] = dwt(x(:,i), base);
    end
    
    for i = 1:N*2
        [dwtLev2(:,i), dwtLev2(:,2*N+i)] = dwt(dwtLev1(:,i), base);
    end
    
    for i = 1:N*4
        [dwtLev3(:,i), dwtLev3(:,4*N+i)] = dwt(dwtLev2(:,i), base);
    end
    
    for i = 1:N*8
        [dwtLev4(:,i), dwtLev4(:,8*N+i)] = dwt(dwtLev3(:,i), base);
    end
    
    for i = 1:N*16
        [dwtLev5(:,i), dwtLev5(:,16*N+i)] = dwt(dwtLev4(:,i), base);
    end
    
    for i = 1:N*32
        [dwtLev6(:,i), dwtLev6(:,32*N+i)] = dwt(dwtLev5(:,i), base);
    end
    
    PCA_feat = zeros(features*2^dwtDepth, N);
    for i = 1:N
        for j = 1:2^dwtDepth
            PCA_feat(features*(j-1)+1,i) = var(dwtLev6(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+2,i) = skewness(dwtLev6(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+3,i) = kurtosis(dwtLev6(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+4,i) = sqrt(sum((dwtLev6(:,N*(j-1)+i)).*dwtLev6(:,N*(j-1)+i))/length(dwtLev6(:,1)));
            PCA_feat(features*(j-1)+5,i) = mean(dwtLev6(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+6,i) = (max(dwtLev6(:,N*(j-1)+i))-min(dwtLev6(:,N*(j-1)+i))) / sqrt(mean([dwtLev6(:,N*(j-1)+i)].^2)); %peak 2 rms. aka crest factor
        end
    end
end

x_ps = x_ps_temp(1:length(x_ps_temp(:,1))/2,:);
x_ps = 20*log10(x_ps);

%normalize and get rid off offset of power spectrum
scaleMean = zeros(length(x_ps(:,1)),2^dwtDepth);
for j =1:2^dwtDepth
    for i = 1:length(x_ps(:,1))
        tempXs = x_ps(i,1+(j-1)*N:j*N);
        scaleMean(i,j) = mean(tempXs);
        x_ps(i,1+(j-1)*N:j*N) = tempXs - mean(tempXs);
    end
end

% normalized std for each spectrum data
for i = 1:length(x_ps(:,1))
    for j = 1:2^dwtDepth
        x_ps(i,(j-1)*N+1:j*N) = x_ps(i,(j-1)*N+1:j*N)./ sqrt(sum(x_ps(i,(j-1)*N+1:j*N).*x_ps(i,(j-1)*N+1:j*N))); 
    end
end


%PCA modeling
x_ps_t = x_ps';
pcsN = 1;
%trans = zeros(N*2^dwtDepth, length(x_ps(:,1)));
trans = zeros(N*2^dwtDepth, pcsN);
Vt = zeros(pcsN*2^dwtDepth, length(x_ps(:,1)));
for i = 1:2^dwtDepth
    [U , S, V] = svd(x_ps_t( 1+(i-1)*N:i*N, : ));
    Vt(1+pcsN*(i-1):pcsN*i, :) = V(:,pcsN)';
    %trans = (pcs'*x_ps)';
    %trans(1+(i-1)*N:i*N,:) = (V'*x_ps( :,1+(i-1)*N:i*N ))';
    trans(1+(i-1)*N:i*N,1:pcsN) = (Vt(1+pcsN*(i-1):pcsN*i, :)*x_ps( :,1+(i-1)*N:i*N ))';
end

% %%-------freq domian process-------
% h = hann(L);
% %h = 1;
% x_ps = zeros(L/2, N);
% x_h = zeros(L, N);
% 
% % average FFT spectrum
% for j = 1:N
%     for i = 1:FrameN
%         x_h(:,j) = x((i-1)*L+1 : L*i, j);
%         x_h(:,j) = x_h(:,j).*h;
%         x_ps_temp = abs(fft(x_h(:,j)));
%         x_ps_temp = x_ps_temp(1:L/2);
%         x_ps(:,j) = x_ps(:,j) + x_ps_temp; 
%     end
% end
% 
% 
% %normalize and get rid off offset of power spectrum
% 
% x_ps = 20*log10(x_ps);
% 
% %get rid of offese
% scaleMean = zeros(L/2,1);
% for i = 1:L/2
%     tempXs = x_ps(i,:);
%     scaleMean(i) = mean(tempXs);
%     x_ps(i,:) = tempXs - mean(tempXs);
% end
% 
% 
% % normalized std for each spectrum data
% for i = 1:L/2
%     x_ps(i,:) = x_ps(i,:)./ sqrt(sum(x_ps(i,:).*x_ps(i,:))); 
% end
% % ------------------end of freq domain process


%PCA Analysis
%[pcs, trans, evs] = princomp(x_ps');
%[pcs, trans, evs] = princomp(x_ps_peak');
%[trans,pcs,evs] = pca1(x_ps_peak);
% trans = real(trans);
% % -----------------Below are SVD transpose in freq
% pcsN = 10;
% [U , S, V] = svd(x_ps');
% %trans = (pcs'*x_ps)';
% trans = (V'*x_ps)';

%test: normalized projection on each pc
% for i = 1:pcsN
% trans(:,i) = trans(:,i)./S(i,i);
% end

for i = 1:num_temp0
    for j = 1:2^dwtDepth
        good(i,j) = trans((i-1)+(j-1)*N+1,1);
    end
end

for i = 1:num_temp1
    for j = 1:2^dwtDepth
        bad(i,j) = trans((i-1)+(j-1)*N+num_temp0+1,1);
    end
end

for i = 1:num_temp2
    for j = 1:2^dwtDepth
        minor(i,j) = trans((i-1)+(j-1)*N+num_temp1+num_temp0+1,1);
    end
end

figure
hold on
for i = 1:num_temp0
    plot(good(i,:))
end

for i = 1:num_temp1
    plot(bad(i,:), 'r')
end

for i = 1:num_temp2
    plot(minor(i,:), 'g')
end



%k mean cluster
%cntPC = 2^dwtDepth;
cntPC = 4 % 1:4 seems more distinguishable
X = zeros(num_temp0 + num_temp1 + num_temp2,cntPC);
X(1:num_temp0,1:cntPC) = good(:,1:1:cntPC);
X(num_temp0+1:num_temp0 + num_temp1,1:cntPC) = bad(:,1:1:cntPC);
X(num_temp0+num_temp1+1:N,1:cntPC) = minor(:,1:1:cntPC);

[idx,c1] = kmeans(X,group)

% Find my cluster centers for group1 2 and 3
% Group1: pass. Group2:fail. 
% Group3: seems fail. Pick files 12 80 96 97 98 109 110 111 112. These files seems like fail 
c2 = zeros(group,cntPC);
c2(1,1:cntPC) = sum(good(:,1:1:cntPC)) / num_temp0;
c2(2,1:cntPC) = sum(bad(:,1:1:cntPC)) / num_temp1;  
c2(3,1:cntPC) = sum(minor(:,1:1:cntPC)) / num_temp2; 
 
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



