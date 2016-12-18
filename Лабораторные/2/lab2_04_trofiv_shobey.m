function [] = lab2_04_trofiv_shobey()

    clc;
    
    n = 10;
    frameSize = 240;
    inputFileName = 'INST.WAV';
    step = 0.00306;
    initinalStepE = 0.001;
    repeats = 30;
    
    compression = zeros(repeats, 1);
    relativeError = zeros(repeats, 1);
    
    for i = 1 : repeats
        
        stepE = initinalStepE * i;
        traceFrame = i;
        
        fprintf('Input file: %s\n', inputFileName);
        fprintf('Step: %f\n', step);
        fprintf('StepE: %f\n', stepE);
        
        [compression(i), relativeError(i)] = process(inputFileName, n, step, stepE, frameSize, traceFrame);
        
        fprintf('\n\n');
        
    end
    
    plot(compression, relativeError);
    
end

function [header, data] = readFile(inputFileName)
    file = fopen(inputFileName, 'r');
    stream = fread(file, 'short');
    fclose(file);
    header = stream(1 : 44);
    data = stream(45 : length(stream));
end

function [] = writeFile(header, data, outputFileName)
    file = fopen(outputFileName, 'w');
    header = header';
    buffer = [header data];
    fwrite(file, buffer, 'short');
    fclose(file);
end

function[result] = entropy(x)

    maxX = max(x);
    minX = min(x);
    result = 0;
    p = zeros(maxX - minX + 1, 1);

    for k = 1 : length(x)
        p(x(k) - minX + 1) = p(x(k) - minX + 1) + 1;
    end
    
    for k = 1 : length(p)
        if p(k) ~= 0
            result = result - (p(k) / length(x) * log2(p(k) / length(x)));
        end
    end
    
end

function [compression, relativeError] = process(inputFileName, m, step, stepE, frameSize, traceFrame)

    [header, data] = readFile(inputFileName);

    frameCount = floor(length(data) / frameSize);
    fprintf('Total frames: %d\n\n', frameCount);

    input = zeros(frameCount * frameSize, 1);
    
    quantizedA = [];
    quantizedE = [];
    
    restoredA = [];
    restoredE = [];
    
    restoredData = [];
    
    squareDifferencefForD = 0;
    squareForD = 0;
    
    p = m + 1;

    for i = 1 : frameCount
    
        if (i == traceFrame) 
            fprintf('Trace for frame %d:\n', i);
        end
   
        currentFrame = data((i - 1) * frameSize + 1 : i * frameSize);
        input((i - 1) * frameSize + 1 : i * frameSize) = currentFrame;
    
        % Autocorrelation method
        c = zeros(p);
        
        for j = 1 : p
            for k = 1 : p
                l = currentFrame(1 : frameSize - abs(j - k));
                r = currentFrame(abs(j - k) + 1 : frameSize);
                c(j, k) = l' * r / frameSize;
            end
        end
        
        r = zeros(1, p);
        for j = 1 : p
            r(j) = c(1, j) / c(1, 1);
        end
        
        if (i == traceFrame)
            fprintf('Yule–Walker coefficients:\n');
            for j = 1 : p
                fprintf('R(%d) = %9.6f\n', j, r(j));
            end
        end
    
        % Levinson-Durbin
        E = zeros(1, p);
        E(1) = r(1);
        a = 0;
        
        for j = 2 : p
            
            if (j == 2)
                currentA = (r(j) - (0 * r(j - 1))) / E(j - 1);
                a = currentA;
            else          
                currentA = (r(j) - sum(a .* r(j - 1 : -1 : 2))) / E(j - 1);
                a(1 : j - 2) = a(1 : j - 2) - currentA * a(j - 2 : -1 : 1);
                a = [a currentA]; %#ok<AGROW>
            end
            
            currentE = E(j - 1) * (1 - a(j - 1) ^ 2);
            E(j) = currentE;
            
        end
       
        if (i == traceFrame)
            fprintf('\nLevinson-Durbin solution:\n');
            for j = 1 : m
                fprintf('a(%d) = %9.6f\n', j, a(j));
            end
            fprintf('E = %9.6f\n\n', E(m));
        end
        
        tempCurrent = zeros(1, frameSize + p - 1);
        tempCurrent(p : length(tempCurrent)) = currentFrame;
        e = zeros(frameSize + p - 1);

        for j = p : frameSize + p - 1
            e(j) = tempCurrent(j) - sum(a .* tempCurrent(j - 1 : -1 : j - p + 1));
        end

        quantizedA = [quantizedA round(a / step)]; %#ok<AGROW>
        temp = round(e / stepE);
        quantizedE = [quantizedE temp(p : frameSize + p - 1)];  %#ok<AGROW>
        
        indexStart = (i - 1) * frameSize + 1;
        indexEnd = (i - 1) * frameSize + frameSize;
        quantizedEBlock = (quantizedE(indexStart : indexEnd));
        
        restoredE = [restoredE stepE * quantizedEBlock]; %#ok<AGROW>
        
        indexStart = (i - 1) * m + 1;
        indexEnd = (i - 1) * m + m;
        quantizedABlock = quantizedA(indexStart : indexEnd);
        
        restoredA = [restoredA step * quantizedABlock]; %#ok<AGROW>
        
        indexStart = (i - 1) * frameSize + 1;
        indexEnd = (i - 1) * frameSize + frameSize;
        tempE = zeros(1, frameSize + p - 1);
        tempE(p : frameSize + p - 1) = restoredE(indexStart : indexEnd);
        
        temp = zeros(1, frameSize + p - 1);
        
        for j = p : frameSize + p - 1
            l = restoredA((i - 1) * m + 1 : (i - 1) * m + m);
            r = temp(j - 1 : -1 : j - m);
            summ = sum(l .* r);
            temp(j) = tempE(j) + summ;
        end

        restoredData = [restoredData temp(p : frameSize + p - 1)]; %#ok<AGROW>

        encoded = temp(p : p + frameSize - 1);
        startIndex = (i - 1) * frameSize + 1;
        endIndex = (i - 1) * frameSize + frameSize;
        original = data(startIndex : endIndex);
        
        for j = 1 : frameSize
            squareForD = squareForD + (encoded(j)) ^ 2;
            squareDifferencefForD = squareDifferencefForD + (encoded(j) - original(j)) ^ 2;
        end

        if (i == traceFrame)
           
            fprintf('Amplitude function\n');
            fprintf('A(w)=sqrt( (1 ');
            
            for alpha = 1 : p - 1
                if (a(alpha) > 0) 
                    fprintf('- %4.3f * cos(%d wT) ', a(alpha), alpha); 
                else
                    fprintf('+ %4.3f * cos(%d wT) ', -a(alpha), alpha);
                end
            end
            
            fprintf(')^2 + ( ');
            
            for alpha = 1 : p - 1
                if (a(alpha) < 0) 
                    fprintf('- %4.3f * sin(%d wT) ', -a(alpha), alpha); 
                else
                    fprintf('+ %4.3f * sin(%d wT) ', a(alpha), alpha); 
                end;
            end;
            
            fprintf(')^2)\n\n');
        
            fprintf('Prediction filter: \ne(n) = x(n)-');
            for j = 1 : m
                fprintf('(%.3f)*x(n-%i)', a(j), j);
                if (j ~= m)
                    fprintf('-');
                end
            end
            fprintf('\n\n');

            fprintf('Synthesis filter: \nxs(n) = e(n)+');
            for j = 1 : m
                fprintf('(%.3f)*x(n-%i)', a(j), j);
                if (j ~= m)
                    fprintf('+');
                end
            end;
            fprintf('\n\n');

            fprintf('Transfer function of prediction filter: \nA(z)=1-((%.3f)*(1/z)+', a(1));
            for j = 2 : m
                fprintf('(%.3f)*(Z^-%i)', a(j),j);
            end;
            fprintf(')\n\n');

            fprintf('Transfer function of synthesis filter: \nH(Z)=1/(1 - ((%.3f)*(1/z)', a(1));
            for j = 2 : p - 1
                fprintf('+(%.3f)*(Z^-%i)', a(j),j);
            end;
            fprintf('))\n\n');
            
        end

    end

    name = strcat('restore', int2str(stepE), '.wav');
    
    writeFile(header, restoredData, name);

    D = squareDifferencefForD / (frameCount * frameSize);
    Drelative = squareDifferencefForD / squareForD;
    fprintf('D: %f\nD relative: %f\n\n', D, Drelative);

    fprintf('Entropy for a: %f\n', entropy(quantizedA));
    fprintf('Average bit count for a: %f\n', entropy(quantizedA) * (length(quantizedA)));

    fprintf('Entropy for e: %f\n', entropy(quantizedE));
    fprintf('Average bit count for a: %f\n\n', entropy(quantizedE) * (length(quantizedE)));

    he = entropy(quantizedE);
    ha = entropy(quantizedA);
    originalSize = 16 * frameSize * frameCount;
    compressedSize = (ha * m * frameCount + he * frameSize * frameCount);
    compressionLevel = originalSize / compressedSize;
    fprintf('Compression level: %f\n\n', compressionLevel);

    relativeError = Drelative;
    compression = compressionLevel;
    
end