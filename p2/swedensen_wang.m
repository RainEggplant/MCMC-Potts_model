function swedensen_wang(K, N_Q, BETA, N, BURN_IN_FACTOR)
disp('Running Swedensen-Wang algorithm:');
tic
n_skip = round(N * BURN_IN_FACTOR); % 统计时跳过最前面样本的个数
x = zeros(K, K, N);
x(:, :, 1) = randi([1, N_Q], K, K);

% 进行迭代
for t = 2:N
    x_t = x(:, :, t - 1);    
    x_cluster = get_cluster(x_t, BETA);
    labels = unique(x_cluster)';
    for l = labels
        x_t(x_cluster == l) = randi([1, N_Q]);
    end

    x(:, :, t) = x_t;  
end

u = calc_u(x(:, :, n_skip:N));
u_mean = mean(u);
toc
disp(['E{u(x)} = ', num2str(u_mean)]);
% ln_z = LN_Z0 - u_mean * BETA;
% disp(['ln Z(T) = ', num2str(ln_z)]);
figure;
histogram(u / (K ^ 2));


    % 采用并查集算法 (Hoshen–Kopelman algorithm) 获取团簇
    function x_cluster = get_cluster(x_t, beta)
        K = length(x_t);
        labels = 1:K ^ 2;
        label_ranks = ones(1, K ^ 2);
        x_cluster = zeros(K, K);
        n_label = 0;

        for row = 1:K
            for col = 1:K      
                % 当左/上方取值与当前值相同时，以 1-e^(-beta) 概率连接
                if col == 1
                    link_left = 0; % 将周期边界留到后面处理
                else
                    link_left = (x_t(row, col) == x_t(row, col - 1) && ...
                        rand() < 1 - exp(-beta));
                end

                if row == 1
                    link_up = 0; % 将周期边界留到后面处理
                else
                    link_up = (x_t(row, col) == x_t(row - 1, col) && ...
                        rand() < 1 - exp(-beta));
                end

                if ~link_left && ~link_up
                    n_label = n_label + 1;
                    x_cluster(row, col) = n_label;
                elseif link_left && ~link_up
                    x_cluster(row, col) = find_label(x_cluster(row, col - 1));
                elseif ~link_left && link_up
                    x_cluster(row, col) = find_label(x_cluster(row - 1, col));
                else
                    x_cluster(row, col) = union_label( ...
                        x_cluster(row, col - 1), x_cluster(row - 1, col));
                end
            end
        end

        % 对周期边界进行处理
        for row = 1:K
            link_left = (x_t(row, 1) == x_t(row, K) && ...
                rand() < 1 - exp(-beta));
            if link_left
                union_label(x_cluster(row, 1), x_cluster(row, K));
            end
        end

        for col = 1:K
            link_up = (x_t(1, col) == x_t(K, col) && ...
                rand() < 1 - exp(-beta));
            if link_up
                union_label(x_cluster(1, col), x_cluster(K, col));
            end
        end

        % 重命名，不过不保证 label 是连续整数
        x_cluster = arrayfun(@(x) find_label(x), x_cluster);


            % 寻找集合的根节点，同时进行路径压缩
            function label = find_label(label_to_find)
                % 寻找根节点
                root = label_to_find;
                while (labels(root) ~= root)
                    root = labels(root);
                end

                % 压缩路径
                while (labels(label_to_find) ~= root)
                    parent = labels(label_to_find);
                    labels(label_to_find) = root;
                    label_to_find = parent;
                end

                label = root;
            end

            % 按树的秩合并两个集合
            function root = union_label(label_1, label_2)
                root = find_label(label_1);
                root_2 = find_label(label_2);
                if root == root_2
                    return;
                end

                if label_ranks(root) > label_ranks(root_2)
                    labels(root_2) = root;
                    label_ranks(root) = ...
                        label_ranks(root) + label_ranks(root_2);
                else
                    labels(root) = root_2;
                    label_ranks(root_2) = ...
                        label_ranks(root) + label_ranks(root_2);
                    root = root_2;
                end
            end

%             function root = union_label(label_1, label_2)
%                 root = find_label(label_1);
%                 root_2 = find_label(label_2);
%                 if root == root_2
%                     return;
%                 end
%                 
%                 labels(root_2) = root;
%             end
    end
end