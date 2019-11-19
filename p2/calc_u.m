% º∆À„ƒ‹¡ø
function u = calc_u(x)
[K, ~] = size(x);
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
