function main ()
clc;
psnry=[];
psnru=[];
psnrv=[];
cz=[];
for i=16:16:128
    [tpsnry tpsnru tpsnrv tcz] = lab3(i);
    fprintf('%f %f %f %f\n',tpsnry,tpsnru,tpsnrv,tcz);
    psnry = [psnry tpsnry];
    psnru = [psnru tpsnru];
    psnrv = [psnrv tpsnrv];
    cz = [cz tcz];    
end
plot(psnry,cz,psnru,cz,psnrv,cz);
end

function [PSNRY,PSNRU,PSNRV,CZ] = lab3(stepAC)
    
    fd = fopen('input.bmp','rb');
    title = fread(fd,54 ,'uint8');
    filename = 'input.bmp';
    %reading into rgb
    global w;
    w = title(19);
    global h;
    h = title(23);

    pixels = fread(fd,[w*3, h],'uint8');
    fclose(fd);

    for i = 1:w
        for j = 1:h
            B(i,j) = pixels((i-1)*3 + 1,j);
            G(i,j) = pixels((i-1)*3 + 2,j);
            R(i,j) = pixels((i-1)*3 + 3,j);
        end
    end

    Y = 0.299*R+0.587*G+0.114*B;
	U = (B-Y)*0.5643 +128;
	V = (R-Y)*0.7132+128;
    
    [resY,PSNRY,CY]=coding(Y,stepAC);
    [resU,PSNRU,CU]=coding(U,stepAC);
    [resV,PSNRV,CV]=coding(V,stepAC);
    CZ=w*h*24/(CY+CU+CV);
    
    fprintf ('C= %f\n',CZ);    
	

    G = resY - 0.714*(resV - 128) - 0.334*(resU - 128);
    R = resY + 1.402*(resV - 128);
    B = resY + 1.772*(resU - 128);
    
    filename = 'restored_lena_';
    ext='.bmp';
    filename = strcat ('restored_lena_',num2str(stepAC),'.bmp');
    file = fopen(filename, 'wb');
    fwrite(file,title,'uint8');
    bytes = zeros(w*3, h);
    for i = 1:w
        for j = 1:h
            bytes((i-1)*3 + 1,j) = B(i,j);
            bytes((i-1)*3 + 2,j) = G(i,j);
            bytes((i-1)*3 + 3,j) = R(i,j);
        end
    end
    fwrite(file,bytes,'uint8');
    fclose(file);
    
end
function[res] = zigzag(X)
    res=[X(1,1),X(1,2),X(2,1),X(3,1),X(2,2),X(1,3),X(1,4),X(2,3),...
       X(3,2),X(4,1),X(5,1),X(4,2),X(3,3),X(2,4),X(1,5),X(1,6),...
       X(2,5),X(3,4),X(4,3),X(5,2),X(6,1),X(7,1),X(6,2),X(5,3),...
       X(4,4),X(3,5),X(2,6),X(1,7),X(1,8),X(2,7),X(3,6),X(4,5),...
       X(5,4),X(6,3),X(7,2),X(8,1),X(8,2),X(7,3),X(6,4),X(5,5),...
       X(4,6),X(3,7),X(2,8),X(3,8),X(4,7),X(5,6),X(6,5),X(7,4),...
       X(8,3),X(8,4),X(7,5),X(6,6),X(5,7),X(4,8),X(5,8),X(6,7),...
       X(7,6),X(8,5),X(8,6),X(7,7),X(6,8),X(7,8),X(8,7),X(8,8)];
end

function[h] = enthropy(x)
a=max(x);
b=min(x);
h=0;
p=zeros(a-b+1);
for k=1:length(x)
    p(x(k)-b+1)=p(x(k)-b+1)+1;
end
for k=1:length(p)
    if p(k)~=0
    h=h-(p(k)/length(x)*log2(p(k)/length(x)));
    end
end
end

function [resY,PSNR,C] = coding(Y,stepAC)
global w h
DC = [];
notAzero = [];
zeroSerLen = [];
notZeroSerLen = [];
stepDC = 8;
resY = [w,h];
    for i=1:w/8
        for j=1:h/8
            block = Y((i-1)*8+1:(i-1)*8+8,(j-1)*8+1:(j-1)*8+8);
            block=dct2(block);
            DC = [DC round((block(1,1)/stepDC))*stepDC];
            block = round(block./stepAC);
        %    DC = [DC block(1,1)*(stepAC/stepDC)];
            AC = zigzag(block);
            zeroCount = 0;
            notZeroCount = 0;
            for g=1:length(AC)
                    if (AC(g)==0)
                        zeroCount = zeroCount+1;
                    else
                        notAzero = [notAzero AC(g)];
                        zeroSerLen = [zeroSerLen zeroCount];
                        zeroCount = 0;
                        notZeroCount = notZeroCount+1;
                    end
            end
            notZeroSerLen = [notZeroSerLen notZeroCount];
            block = block.*stepAC;
            block = idct2(block);
            
            resY((i-1)*8+1:(i-1)*8+8,(j-1)*8+1:(j-1)*8+8) = block;
        end
    end

    
    DC = DC(2:length(DC)) - DC(1:length(DC)-1);
    C = enthropy(DC)*(h*w/64)+enthropy(zeroSerLen)*length(zeroSerLen)+enthropy(notAzero)*length(notAzero)+enthropy(notZeroSerLen)*(h*w/64);
    PSNR = 10*log10((255)^2/sum(sum((Y-resY).*(Y-resY)))*w*h);
 


    SNR=10*log10((255)^2/sum(sum((Y-resY).*(Y-resY)))*w*h);

 end

                       
    
