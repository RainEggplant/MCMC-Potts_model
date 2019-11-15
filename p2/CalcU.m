% º∆À„ƒ‹¡ø
function u = CalcU(x)
[K, ~, T] = size(x);
u = zeros(T, 1);

for t = 1:T
    for row = 1:K
       for col = 1:K
          if x(row, col, t) == x(mod(row, K) + 1, col, t) 
             u(t) = u(t) - 1; 
          end
          if x(row, col, t) == x(row, mod(col, K) + 1, t) 
             u(t) = u(t) - 1; 
          end
       end
    end
end
end
