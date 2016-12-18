function lab1_04_trofiv_shobey()
    clc;
    x = [-26.8932 -56.5759 -18.5382 -10.5765 -4.0511 -0.8030 33.1166 54.0909 -5.9260 1.9893 -2.4746 -1.6147 0.0884 -50.4010 3.2252 -1.6825 1.4515 4.2397 1.3388 1.4655 15.2490 -3.0892 0.4293 -1.7859 -2.4400 0.3259 -1.1790 -12.0110 -2.7748 1.3800 5.4677 0.1175 1.9084 0.9616 6.6647 -1.3044 -1.8629 -1.2685 -0.2629 -0.5833 -0.1081 -2.4758 0.5333 0.7980 -0.4441 0.1581 0.7309 -0.2500 1.9203 -0.4097 -0.1532 -0.8937 0.6671 0.1928 -0.3112 0.7015 0.6157 -0.4927 -1.2141 -0.3262 -0.2809 0.7488 0.1304];
    step = 6;
    
    xmax = max(x);
    xmin = min(x);
    r  = max(abs(xmax), abs(xmin));
    
    fprintf('Входная последовательность:\n');
    fprintf('%.4f ', x);
    
    fprintf('\nКвантованная последовательность с переменной скоростью с шагом %d:\n', step);
    y = round(x / step);
    fprintf('%d ', y);
    
    fprintf('\nВосстановленная последовательность:\n');
    z = y * step;
    fprintf('%.4f ', z);
    
    d = 0;
    e = 0;
    for i = 1 : length(x)
        d = d + (abs(x(i)) - abs(z(i)))^2;
        e = e + x(i)^2;
    end
    d = d / length(x);
    e = e / length(x);
    fprintf('\nОтносительная среднеквадратичная ошибка:\n');
    fprintf('%.4f', d / e);
    
    fprintf('\nВероятности аппроксимирующих значений:\n');
    w = unique(y);
    p = zeros(1, length(w));
    for i = 1 : length(y)
        for j = 1 : length(w)
            if (y(i) ==  w(j))
                p(j) = p(j) + 1;
                break;
            end
        end
    end
    for i = 1 : length(p)
        p(i) = p(i) / length(z);
        fprintf('%3d: %.5f\n', w(i), p(i));
    end
    
    fprintf('Энтропия на выходе квантователя:\n');
    h = 0;
    for i = 1 : length(p)
        h = h - p(i) * log2(p(i));
    end
    fprintf('%.5f\n', h);
    
    fprintf('Среднее число битов для хранения закодированной последовательности:\n');
    fprintf('%d\n', round(h * length(y)));
    
    fprintf('Квантованная c помощью 4 квантов последовательность:\n');
    q = 4;
    step = 2 * r / q;
    quantum = zeros(1, q);
    for i = 1 : q
        quantum(i) = -r + (i - 1) * step + step / 2; 
    end
    for i = 1 : length(x)
        for j = 1 : q
            diff = abs(x(i) - quantum(q - j + 1));
            if (diff <= step / 2)
                z(i) = quantum(q - j + 1);
                y(i) = q - j;
                break;
            end
        end
    end
    fprintf('%d ', y);
    fprintf('\nВосстановленная последовательность:\n');
    fprintf('%.3f ', z);
    fprintf('\nЧисло битов для хранения закодированной c помощью 4 квантов последовательности:\n');
    fprintf('%d', round(length(x) * log2(q)));
    d = 0;
    for i = 1 : length(x)
        d = d + (abs(x(i)) - abs(z(i)))^2;
    end
    d = d / length(x);
    fprintf('\nОтносительная среднеквадратичная ошибка:\n');
    fprintf('%.4f\n', d / e);
    
    fprintf('Квантованная c помощью 8 квантов последовательность:\n');
    q = 8;
    step = 2 * r / q;
    quantum = zeros(1, q);
    for i = 1 : q
        quantum(i) = -r + (i - 1) * step + step / 2; 
    end
    for i = 1 : length(x)
        for j = 1 : q
            diff = abs(x(i) - quantum(q - j + 1));
            if (diff <= step / 2)
                z(i) = quantum(q - j + 1);
                y(i) = q - j;
                break;
            end
        end
    end
    fprintf('%d ', y);
    fprintf('\nВосстановленная последовательность:\n');
    fprintf('%.3f ', z);
    fprintf('\nЧисло битов для хранения закодированной c помощью 4 квантов последовательности:\n');
    fprintf('%d', round(length(x) * log2(q)));
    d = 0;
    for i = 1 : length(x)
        d = d + (abs(x(i)) - abs(z(i)))^2;
    end
    d = d / length(x);
    fprintf('\nОтносительная среднеквадратичная ошибка:\n');
    fprintf('%.4f\n', d / e);
end