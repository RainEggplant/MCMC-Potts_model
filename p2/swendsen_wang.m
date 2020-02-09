function [u, x] = swendsen_wang(K, N_Q, BETA, N)
x = randi([1, N_Q], K, K); % ��һ���������ȡֵ
u = zeros(N, 1); % ����
u(1) = calc_u(x);

% ���е���
for t = 2:N
    x_cluster = get_cluster(x, BETA);
    labels = unique(x_cluster)';
    for l = labels
        x(x_cluster == l) = randi([1, N_Q]);
    end
    u(t) = calc_u(x);
end

% ln_z = LN_Z0 - u_mean * BETA;
% disp(['ln Z(T) = ', num2str(ln_z)]);
% figure;
% histogram(u / (K ^ 2));


    % ���ò��鼯�㷨 (Hoshen�CKopelman algorithm) ��ȡ�Ŵ�
    function x_cluster = get_cluster(x_t, beta)
        K = length(x_t);
        labels = 1:K ^ 2;
        label_ranks = ones(1, K ^ 2);
        x_cluster = zeros(K, K);
        n_label = 0;

        for row = 1:K
            for col = 1:K      
                % ����/�Ϸ�ȡֵ�뵱ǰֵ��ͬʱ���� 1-e^(-beta) ��������
                if col == 1
                    link_left = 0; % �����ڱ߽��������洦��
                else
                    link_left = (x_t(row, col) == x_t(row, col - 1) && ...
                        rand() < 1 - exp(-beta));
                end

                if row == 1
                    link_up = 0; % �����ڱ߽��������洦��
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

        % �����ڱ߽���д���
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

        % ����������������֤ label ����������
        x_cluster = arrayfun(@(x) find_label(x), x_cluster);


            % Ѱ�Ҽ��ϵĸ��ڵ㣬ͬʱ����·��ѹ��
            function label = find_label(label_to_find)
                % Ѱ�Ҹ��ڵ�
                root = label_to_find;
                while (labels(root) ~= root)
                    root = labels(root);
                end

                % ѹ��·��
                while (labels(label_to_find) ~= root)
                    parent = labels(label_to_find);
                    labels(label_to_find) = root;
                    label_to_find = parent;
                end

                label = root;
            end

            % �������Ⱥϲ���������
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

    % ��������
    function u = calc_u(x)
        u = 0;

        for row = 1:K
            for col = 1:K
                if x(row, col) == x(row, mod(col, K) + 1) 
                    u = u - 1;
                end
                if x(row, col) == x(mod(row, K) + 1, col) 
                    u = u - 1;
                end
            end
        end
    end
end
