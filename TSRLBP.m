function TSRLBPfeatures = TSRLBP(Image,radius,neighbors,mappMode,histMode,Coding,combMode,table,newMax,block_size)

if size(Image,3) == 3
    Image = rgb2gray(Image);
    Image = double(Image);
else
    Image = double(Image);
end

% mapmode normal,u2,ri and riu2
samples = neighbors;
% [table,newMax] = getmapping(samples,mappMode);

% Coding: 64 different Codings
tableCode0 = fliplr(de2bi(Coding));
tableCode00 = zeros(1,6-length(tableCode0));
tableCode1 = fliplr([tableCode00,tableCode0]);

d_image = double(Image);
spoints = zeros(neighbors,2);

% Angle step.
a = 2*pi/neighbors;
for i = 1:neighbors
    spoints(i,1) = -radius*sin((i-1)*a);
    spoints(i,2) = radius*cos((i-1)*a);
end

% Determine the dimensions of the input image.
[ysize xsize] = size(Image);
miny = min(spoints(:,1));
maxy = max(spoints(:,1));
minx = min(spoints(:,2));
maxx = max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey = ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex = ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block
origy = 1-floor(min(miny,0));
origx = 1-floor(min(minx,0));

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.
Center = Image(origy:origy+dy,origx:origx+dx);
bins = 2^neighbors;


Intensitys = zeros(dy+1,dx+1,neighbors);
mu0 = zeros(dy+1,dx+1,neighbors);

for i = 1:neighbors
    y = spoints(i,1)+origy;
    x = spoints(i,2)+origx;
    % Calculate floors, ceils and rounds for the x and y.
    fy = floor(y); cy = ceil(y); ry = round(y);
    fx = floor(x); cx = ceil(x); rx = round(x);
    % Check if interpolation is needed.
    if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
        % Interpolation is not needed, use original datatypes
        N = Image(ry:ry+dy,rx:rx+dx);
    else
        % Interpolation needed, use double type images 
        ty = y-fy;
        tx = x-fx;

        % Calculate the interpolation weights.
        w1 = roundn((1-tx)*(1-ty),-6);
        w2 = roundn(tx*(1-ty),-6);
        w3 = roundn((1-tx)*ty,-6) ;
        % w4 = roundn(tx * ty,-6) ;
        w4 = roundn(1-w1-w2-w3,-6);

        % Compute interpolated pixel values
        N = w1*d_image(fy:fy+dy,fx:fx+dx)+w2*d_image(fy:fy+dy,cx:cx+dx) + ...
            w3*d_image(cy:cy+dy,fx:fx+dx)+w4*d_image(cy:cy+dy,cx:cx+dx);
        N = roundn(N,-4);
    end
    mu0(:,:,i) = abs(N-Center);
    Intensitys(:,:,i) = N;
end


indN = [1:neighbors,1:neighbors];
for jk = 1:neighbors
    pn(:,:,:,jk) = cat(3,Center,Intensitys(:,:,indN(jk+combMode)),Intensitys(:,:,jk));
end

for jk = 1:neighbors
    code1 = double(((pn(:,:,1,jk)>=pn(:,:,2,jk)) .* (pn(:,:,1,jk)>=pn(:,:,3,jk)) .* (pn(:,:,2,jk)>=pn(:,:,3,jk))));
    c2 = double((code1==0));
    code2 = c2 .* double(((pn(:,:,1,jk)>=pn(:,:,2,jk)) .* (pn(:,:,1,jk)>=pn(:,:,3,jk)) .* (pn(:,:,3,jk)>=pn(:,:,2,jk))));
    c3 = double((code1==0).*(code2==0));
    code3 = c3 .* double(((pn(:,:,2,jk)>=pn(:,:,1,jk)) .* (pn(:,:,2,jk)>=pn(:,:,3,jk)) .* (pn(:,:,1,jk)>=pn(:,:,3,jk))));
    c4 = double((code1==0).*(code2==0).*(code3==0));
    code4 = c4 .* double(((pn(:,:,2,jk)>=pn(:,:,1,jk)) .* (pn(:,:,2,jk)>=pn(:,:,3,jk)) .* (pn(:,:,3,jk)>=pn(:,:,1,jk))));
    c5 = double((code1==0).*(code2==0).*(code3==0).*(code4==0));
    code5 = c5 .* double(((pn(:,:,3,jk)>=pn(:,:,1,jk)) .* (pn(:,:,3,jk)>=pn(:,:,2,jk)) .* (pn(:,:,2,jk)>=pn(:,:,1,jk))));
    c6 = double((code1==0).*(code2==0).*(code3==0).*(code4==0).*(code5==0));
    code6 = c6 .* double(((pn(:,:,3,jk)>=pn(:,:,1,jk)) .* (pn(:,:,3,jk)>=pn(:,:,2,jk)) .* (pn(:,:,1,jk)>=pn(:,:,2,jk))));
    code_p(:,:,jk) = tableCode1(1)*code1+tableCode1(2)*code2+tableCode1(3)*code3+...
    tableCode1(4)*code4+tableCode1(5)*code5+tableCode1(6)*code6;
end
    
TSRLBP_S0 = 0;
for jk = 1:neighbors
    TSRLBP_S0 = TSRLBP_S0+2^(jk-1)*code_p(:,:,jk);
end


mu = zeros(dy+1,dx+1);
for i = 1:neighbors
    mu = mu+abs(mu0(:,:,i));
end

mu = mu/neighbors;

sig = zeros(dy+1,dx+1);
for i = 1:neighbors
    sig = sig+(mu0(:,:,i)-mu).^2;
end
sig = sig/neighbors;

%Calculate t_mu and t_sig
[A B] = size(Center);
mu_t = zeros(A,B);
sig_t = zeros(A,B);
mu_I = zeros(A,B);
t1 = fix(A/block_size);
t2 = fix(B/block_size);

for i = 1:block_size
    for j = 1:block_size
        if i<block_size && j<block_size
            mat1 = (mu(1+(i-1)*t1:i*t1,1+(j-1)*t2:j*t2));
            mu_t(1+(i-1)*t1:i*t1,1+(j-1)*t2:j*t2) = mean(mat1(:));
            mat2 = (sig(1+(i-1)*t1:i*t1,1+(j-1)*t2:j*t2));
            sig_t(1+(i-1)*t1:i*t1,1+(j-1)*t2:j*t2) = mean(mat2(:));
            mat3 = (Center(1+(i-1)*t1:i*t1,1+(j-1)*t2:j*t2));
            mu_I(1+(i-1)*t1:i*t1,1+(j-1)*t2:j*t2) = mean(mat3(:));
        end
        if i<block_size & j==block_size
            mat1 = (mu(1+(i-1)*t1:i*t1,1+(j-1)*t2:B));
            mu_t(1+(i-1)*t1:i*t1,1+(j-1)*t2:B) = mean(mat1(:));
            mat2 = (sig(1+(i-1)*t1:i*t1,1+(j-1)*t2:B));
            sig_t(1+(i-1)*t1:i*t1,1+(j-1)*t2:B) = mean(mat2(:));
            mat3 = (Center(1+(i-1)*t1:i*t1,1+(j-1)*t2:B));
            mu_I(1+(i-1)*t1:i*t1,1+(j-1)*t2:B) = mean(mat3(:));
        end
        if i==block_size & j<block_size
            mat1 = (mu(1+(i-1)*t1:A,1+(j-1)*t2:j*t2));
            mu_t(1+(i-1)*t1:A,1+(j-1)*t2:j*t2) = mean(mat1(:));
            mat2 = (sig(1+(i-1)*t1:A,1+(j-1)*t2:j*t2));
            sig_t(1+(i-1)*t1:A,1+(j-1)*t2:j*t2) = mean(mat2(:));
            mat3 = (Center(1+(i-1)*t1:A,1+(j-1)*t2:j*t2));
            mu_I(1+(i-1)*t1:A,1+(j-1)*t2:j*t2) = mean(mat3(:));
        end
        if i==block_size & j==block_size
            mat1 = (mu(1+(i-1)*t1:A,1+(j-1)*t2:B));
            mu_t(1+(i-1)*t1:A,1+(j-1)*t2:B) = mean(mat1(:));
            mat2 = (sig(1+(i-1)*t1:A,1+(j-1)*t2:B));
            sig_t(1+(i-1)*t1:A,1+(j-1)*t2:B) = mean(mat2(:));
            mat3 = (Center(1+(i-1)*t1:A,1+(j-1)*t2:B));
            mu_I(1+(i-1)*t1:A,1+(j-1)*t2:B) = mean(mat3(:));
        end
    end
end

TSRLBP_mu = mu>=mu_t;
TSRLBP_sig = sig>=sig_t;
TSRLBP_I = Center>mu_I;

%Apply mapping if it is defined
TSRLBP_S = zeros(size(TSRLBP_S0));
bins = newMax;
for i = 1:size(TSRLBP_S0,1)
    for j = 1:size(TSRLBP_S0,2)
        TSRLBP_S(i,j) = table(TSRLBP_S0(i,j)+1);
    end
end

TSRLBPfeatures0 = 2^3*TSRLBP_S+2^2*TSRLBP_mu+2^1*TSRLBP_sig+2^0*TSRLBP_I;

bins = newMax*4*2;

if (strcmp(histMode,'h') || strcmp(histMode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    TSRLBPfeatures = hist(TSRLBPfeatures0(:),0:(bins-1));
    if (strcmp(histMode,'nh'))
        TSRLBPfeatures = TSRLBPfeatures/sum(TSRLBPfeatures);
    end
end
