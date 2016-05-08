%MVDR
%author supreet 4/23/2016
%% %% Simulation data
simulatedata = 1; %input 0 = false or 1 = true
timeORfreq = 0; %input 1 = time 0 = freq;
barORcapon = 1; % input 1 = barlett 0 = Capon (MVDR)
sen3 = [1,0.24]; % sensor to localize
Cx = 100;
Cy = 100;
length = 1.22;
bredth = 1.22;
cx = linspace(0, length, Cx);
cy = linspace(0, bredth, Cy);
[CLx, CLy] = meshgrid(cx,cy);
Cxy = [reshape(CLx, Cx*Cy, 1) reshape(CLy, Cx*Cy, 1)];
points = Cxy;
%%
if(simulatedata == 0)
    a = 0.16;
    b = 0.86;
    Tx = [a, b];
    distancee = zeros(size(points,1),1);
    for v = 1:size(points,1)
        distancee(v) = sqrt((Tx(1,1)- points(v,1))^2 + (Tx(1,2) - points(v,2))^2);
    end
    Fs = 10e6;
    Q = 1e4;
    d = distancee;
    [x1, ~, ~] = simlamb(d, Fs, Q, [1]);
    save sen2 x1;
    %%
    a = 0.38;
    b = 0.27;
    Tx = [a, b];
    distancee = zeros(size(points,1),1);
    for v = 1:size(points,1)
        distancee(v) = sqrt((Tx(1,1)- points(v,1))^2 + (Tx(1,2) - points(v,2))^2);
    end
    Fs = 10e6;
    Q = 1e4;
    d = distancee;
    [x2, ~, ~] = simlamb(d, Fs, Q, [1]);
    save sen4 x2;
    
    %%
    a = 0.53;
    b = 0.14;
    Tx = [a, b];
    distancee = zeros(size(points,1),1);
    for v = 1:size(points,1)
        distancee(v) = sqrt((Tx(1,1)- points(v,1))^2 + (Tx(1,2) - points(v,2))^2);
    end
    Fs = 10e6;
    Q = 1e4;
    d = distancee;
    [x3, ~, ~] = simlamb(d, Fs, Q, [1]);
    save sen6 x3;
    %%
    a = 0.61;
    b = 0.97;
    Tx = [a, b];
    distancee = zeros(size(points,1),1);
    for v = 1:size(points,1)
        distancee(v) = sqrt((Tx(1,1)- points(v,1))^2 + (Tx(1,2) - points(v,2))^2);
    end
    Fs = 10e6;
    Q = 1e4;
    d = distancee;
    [x4, ~, ~] = simlamb(d, Fs, Q, [1]);
    save sen8 x4;
    %%
    a = 0.71;
    b = 0.2;
    Tx = [a, b];
    distancee = zeros(size(points,1),1);
    for v = 1:size(points,1)
        distancee(v) = sqrt((Tx(1,1)- points(v,1))^2 + (Tx(1,2) - points(v,2))^2);
    end
    Fs = 10e6;
    Q = 1e4;
    d = distancee;
    [x5, ~, ~] = simlamb(d, Fs, Q, [1]);
    save sen10 x5;
    
    %%
    a = 0.8;
    b = 0.94;
    Tx = [a, b];
    distancee = zeros(size(points,1),1);
    for v = 1:size(points,1)
        distancee(v) = sqrt((Tx(1,1)- points(v,1))^2 + (Tx(1,2) - points(v,2))^2);
    end
    Fs = 10e6;
    Q = 1e4;
    d = distancee;
    [x6, ~, ~] = simlamb(d, Fs, Q, [1]);
    save sen12 x6;
    %%
    a = 0.2;
    b = 0.5;
    Tx = [a, b];
    distancee = zeros(size(points,1),1);
    for v = 1:size(points,1)
        distancee(v) = sqrt((Tx(1,1)- points(v,1))^2 + (Tx(1,2) - points(v,2))^2);
    end
    Fs = 10e6;
    Q = 1e4;
    d = distancee;
    [x7, ~, ~] = simlamb(d, Fs, Q, [1]);
    save senRand x7;
end
%% get simulated observations

for p = 1:size(pointsSen3,1)
    di(p) = sqrt((sen3(1,1)- pointsSen3(p,1))^2 + (sen3(1,2) - pointsSen3(p,2))^2);
end
Fs = 10e6;
Q = 1e4;
d = di;
[obs, ~, ~] = simlamb(d, Fs, Q, [1]);

%% add noise to the simulated observations
variancee = 0.000000004;

for j = 1:size(obs,2)
    obsNoise(:,j) = obs(:,j)+sqrt(variancee)*rand(size(obs,1),1);
end

% for j = 1:size(obs,2)
%     obsNoise(:,j) = obs(:,j)+rand(1,1)*rand(size(obs,1),1);
% end
%% Form observation vectors
if (timeORfreq ==0)
    fn = 20:50;
    obsNoiseF = fft(obsNoise);
    F = obsNoiseF(fn,:);
    %F = fft(obs(fn,1:6));
    F = F./norm(F);
    % construct the spectrum
    R = zeros(size(fn,2),size(fn,2));
    for i = 1:size(F,2)
        R =R + (F(:,i)*(ctranspose(F(:,i))));
    end
    %R = F*ctranspose(F);
    R = R./size(F,2);
    % iterate across all locations
    tm = 0;
    L = 10000;
    pxl1 = zeros(1,L);
    fprintf(repmat(' ',1,41));
    G = eye(size(R));
    Rl = R + 0.0004*G;
    for n = 1:L
        fprintf([ repmat('\b', 1, 41) '%08i / %08i [Time left: %s]'], n, L, datestr(tm/24/3600*(L-n+1), 'HH:MM:SS')); ts = tic;
        vec = zeros(5,1);
        XX_sen2  =  fft(x1(:,n)./(norm(x1(:,n))));
        XX_sen4  =  fft(x2(:,n)./(norm(x2(:,n))));
        XX_sen6  =  fft(x3(:,n)./(norm(x3(:,n))));
        XX_sen8  =  fft(x4(:,n)./(norm(x4(:,n))));
        XX_sen10 =  fft(x5(:,n)./(norm(x5(:,n))));
        XX_sen12 =  fft(x6(:,n)./(norm(x6(:,n))));
%         XX_senrand =  fft(x7(:,n)./(norm(x7(:,n))));
        h =  horzcat(XX_sen2(fn),XX_sen4(fn),XX_sen6(fn),XX_sen8(fn),XX_sen10(fn),XX_sen12(fn));
        h = h./norm(h);
        if (barORcapon == 0)
        pxl1(n) = 1./trace(ctranspose(h)*(R\h)); %MVDR freq
       % pxl(n) = (ctranspose(h)*inv(r))./trace(ctranspose(h)*(R\h)))
        end
        if (barORcapon == 1)
            pxl1(n) = trace(ctranspose(h)*Rl*h); % Barlett freq
        end
        tm = (toc(ts) + tm)/min([n 2]);
    end
    pxl1 = pxl1 - min(pxl1);                   % Normalize -- make minimum value zero
    pxl1 = pxl1 / max(pxl1);                   % Normalize -- make maximum value one
    pxl10 = reshape(pxl1, Cx, Cy);  % Shape into grid
end

%% time domain
if (timeORfreq == 1)
    F = obsNoise(1:2000,:);
    F = F./norm(F);
    Rx = zeros(size(F,1),size(F,1));
    for i = 1:size(F,2)
        Rx =Rx + (F(:,i)*(F(:,i)'));
    end
    L = 10000; tm = 0;
    G = eye(size(Rx));
    Rxl = Rx + 0.002*G;  
    for n = 1:L
        fprintf([ repmat('\b', 1, 41) '%08i / %08i [Time left: %s]'], n, L, datestr(tm/24/3600*(L-n+1), 'HH:MM:SS')); ts = tic;
        vec = zeros(5,1);
        XX_sen2  =  (x1(1:2000,n))./norm(fft(x1(1:2000,n)));
        XX_sen4  =  (x2(1:2000,n))./norm(fft(x2(1:2000,n)));
        XX_sen6  =  (x3(1:2000,n))./norm(fft(x3(1:2000,n)));
        XX_sen8  =  (x4(1:2000,n))./norm(fft(x4(1:2000,n)));
        XX_sen10 =  (x5(1:2000,n))./norm(fft(x5(1:2000,n)));
        XX_sen12 =  (x6(1:2000,n))./norm(fft(x6(1:2000,n)));
        h =  horzcat(XX_sen2,XX_sen4,XX_sen6,XX_sen8,XX_sen10,XX_sen12);
        pxl1(n) = 1./abs(trace(ctranspose(h)*(inv(Rxl))*h)); %MVDR time
        tm = (toc(ts) + tm)/min([n 2]);
    end
    pxl1 = pxl1 - min(pxl1);                   % Normalize -- make minimum value zero
    pxl1 = pxl1 / max(pxl1);                   % Normalize -- make maximum value one
    pxl10 = reshape(pxl1, Cx, Cy);  % Shape into grid
    
end
    
%% plot results
imagesc(cx, cy, abs(pxl10));
hold on;
axis([0 1.22 0 1.22]); axis xy; axis square;
xlabel('Plate width [m]'); ylabel('Plate length [m]');
plot(pointsSen3(:,1),pointsSen3(:,2), 'rx', 'linewidth', 1.5, 'markersize',7);
plot(sen3(:,1), sen3(:,2), 'ks', 'linewidth',2, 'markersize',5);
%plot(sen5(:,1), sen5(:,2), 'ks', 'linewidth',2, 'markersize',5);
hold off;
%legend('Sensor location', 'Source Location', 'Location', 'SouthOutside'); legend BOXOFF;
if(barORcapon == 0)
title('MVDR ALGORITHM or Capon Processor');
end

if (barORcapon == 1)
title('Barlett Processor');
end
    