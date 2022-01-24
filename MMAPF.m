% Min-Max Average Pooling Filter for Impulse Noise Removal
% Takes nImg <Noisy Image> as input

function OutImg = MMAPF(nImg)

% calculation of size of Image
[row, col] = size(nImg);

% Variable to Count number of Noisy Pixels
noisy_pixel_count = 0;
for i=1:row
    for j=1:col
        if (nImg(i,j) == 0) || (nImg(i,j) == 255)
            noisy_pixel_count = noisy_pixel_count + 1;
            b_f(i,j) = 0;
        else
            b_f(i,j) = 1;
        end
    end
end

% Noise Density as ratio of suspected corrupted pixels by total pixels
noise_density = noisy_pixel_count ./ (row * col);

% Initial Noise Density-Dependent Fcuntion for Noise Removal below a
% certain threshold

if (noise_density < 0.45)
    Img_IEHCLND = IEHCLND(nImg, noise_density);
else
    Img_IEHCLND = nImg;
end

I_1 = double(Img_IEHCLND);
I_2 = double(Img_IEHCLND);

% Second stage
[I_1, I_2] = CCMP(I_1, I_2, 'Max'); 
[I_1, I_2] = CCMP(I_1, I_2, 'Min'); 
[I_1, I_2] = CCMP(I_1, I_2, 'Min'); 
[I_1, I_2] = CCMP(I_1, I_2, 'Max'); 

% Final Step for Output Image
OutImg = RnS(I_1, I_2, b_f);

end

function out = IEHCLND(nImg,N_d)
    [row, col] = size(nImg);
    imgZP = padarray(nImg,[1 1]);
    a = floor(N_d./0.1);
    for i = 2 : row +1
        for j = 2 : col + 1
            if (imgZP(i,j) ==0) || (imgZP(i,j) == 255)
                tmp = imgZP(i-1:i+1,j-1:j+1);
                tmp(tmp==0) = [];tmp(tmp==255) = []; %removing noisy pixe
                if length(tmp) > a 
                    pixel = median(tmp);
                else
                    pixel = imgZP(i,j);
                end
            else
                pixel = imgZP(i,j);
            end
            output(i-1,j-1) = uint8(pixel);
        end
    end
    out = output;
end

function [O_1, O_2] = CCMP(I_1, I_2, str)

[row,col] = size(I_1);

I_1 = padarray(I_1, [1 1]);
I_2 = padarray(I_2, [1 1]);

for i = 1 : row + 2
    for j = 1 : col + 2
        if I_1(i,j)==0||I_1(i,j)==255
            b_g(i,j) = 0;
        else
            b_g(i,j) = 1;
        end
    end
end

for i = 2 : row + 1
    for j = 2 : col + 1
        if b_g(i,j) == 0
            temp = nonzeros(I_1(i-1:i+1,j-1:j+1).*b_g(i-1:i+1,j-1:j+1));
            if ~isempty(temp)
                if str == 'Max'
                    I_1(i,j) = max(temp);
                    I_2(i,j) = min(temp);
                elseif str == 'Min'
                    I_1(i,j) = min(temp);
                    I_2(i,j) = max(temp);
                end
            end
        end
    end
end

O_1 = I_1(2:row+1,2:col+1);
O_2 = I_2(2:row+1,2:col+1);

end

function out = RnS(I_1, I_2, b_f)
    oImg = uint8(I_1.* 0.5 + I_2.* 0.5);
    [row, col] = size(oImg);
    for i = 2 : row - 1
        for j = 2 : col - 1
            if b_f(i,j)==0
                oImg(i,j) = round(mean(nonzeros(oImg(i-1:i+1,j-1:j+1)),'all'));
            end
        end
    end
    out = oImg;
end