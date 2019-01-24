function [sol, i, res] = steepest_descent(f, g, H, x0, iterations, search, error)
    res = 0;
    x = x0;
    x_old = x;
    if(strcmp(search, 'fibonacci'))
        fibo = fibonacci(50);
    end
    loss = [];
    
    for i = 1:iterations
        loss = [loss f(x)];
        d_a = -g(x);
        f_a = @(a) f(x + a.*d_a);
%         g_a = g(x);
        alpha_max = a_max_calculate(f, -g(x), x);
        if(strcmp(search, 'lagrange'))
            alpha = lagrange_search(f_a, 0, alpha_max);
          elseif(strcmp(search, 'golden') )
            [alpha, it] = golden_section(f_a, 0, alpha_max, 1e-4);
          elseif(strcmp(search, 'fibonacci'))
            [alpha, it] = fibonacci_search(f_a, 0, alpha_max, fibo, 1e-4);
          elseif(strcmp(search, 'dichotomous'))
            [alpha, it] = dichotomous_search(f_a, 0, alpha_max, 1e-4);
          elseif(strcmp(search, 'three_point'))
            [alpha, it] = threePointInterval(f_a, 0, alpha_max, 1e-4);
          elseif(strcmp(search, 'newton_line'))
            g_a = @(a) d_a'*g(x + a.*d_a);
            H_a = @(a) d_a'*H(x + a.*d_a)*d_a; 
            [alpha, it] = newton_line(f_a, g_a, H_a, 0, alpha_max, 1e-4);
          end
%         alpha = rand();
%        display(alpha);
        x = x - alpha*g(x);
        if(norm(x - x_old) < error)
            res = 1;
%             break;
        end
        x_old = x;
    end
    sol = x;
    i = i+1;
end