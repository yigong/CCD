load('./pos_1.mat');
N = 30;
indices = int16(linspace(1, length(alpha), N));

for i = indices
    plot(ridgeWidth{i,2}, 'b')
    hold on
end

figure(1)
load('./pos_5.mat')
indices = int16(linspace(1, length(alpha), N));

for i = indices
    plot(ridgeWidth{i,2}, 'r')
    hold on
end

figure(1)
load('./pos_9.mat')
indices = int16(linspace(1, length(alpha), N));

for i = indices
    plot(ridgeWidth{i,2}, 'g')
    hold on
end

figure(1)
load('./pos_13.mat')
indices = int16(linspace(1, length(alpha), N));

for i = indices
    plot(ridgeWidth{i,2}, 'c')
    hold on
end