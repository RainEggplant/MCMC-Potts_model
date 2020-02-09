function plot_sample(x, N_Q)
x = x / N_Q;
figure;

for k = 1:5
    subplot(2, 3, k);
    imshow(x(:, :, k));
    title(['T_', num2str(k), ' µäÐÍÑù±¾']);
end
end

