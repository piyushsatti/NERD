function OutImg = IBLF(nImg)

    % Intensity Bound Limit Filter for Impulse Noise Removal
    % Takes nImg <Noisy Image> as input
    
    % calculation of size of Image
    [row, col] = size(nImg);
    
    % flags the noisy elements with 0
    flag = ones([row, col]);
    for i=1:row
        for j=1:col
            if (nImg(i,j) == 0) || (nImg(i,j) == 255)
                flag(i,j) = 0;
            end
        end
    end
    
    % Padding of Image
    pad = 4;
    nImg = double(padarray(nImg, [pad pad]));
    flag = padarray(flag, [pad pad]);
    
    % Declaring Lower and Upper Bound Image Variables
    I_l = zeros([row, col]);
    I_u = zeros([row, col]);
    
    for i = 1 + pad : row + pad
        for j = 1 + pad : col + pad
            if flag(i,j) == 0                
                
                % Default Bound Image Contruction Step
                temp_cross = [reshape(nImg(i-1:i+1,j).*flag(i-1:i+1,j), 1, []), nImg(i,j-1).*flag(i,j-1), nImg(i,j+1).*flag(i,j+1)];
                temp_cross(temp_cross == 0) = [];
                if ~isempty(temp_cross)
                    I_l(i-pad,j-pad) = min(temp_cross);
                    I_u(i-pad,j-pad) = max(temp_cross);
                    continue;
                end
                
                % Adaptive Bound Image Contruction in case default fails
                [I_l(i-pad,j-pad), I_u(i-pad,j-pad)] = ABIC(pad, nImg(i-pad:i+pad,j-pad:j+pad), flag(i-pad:i+pad,j-pad:j+pad));
                
            else
                I_l(i-pad,j-pad) = nImg(i,j);
                I_u(i-pad,j-pad) = nImg(i,j);
            end
        end
    end
    
    % Bound Limit Step
    [I_l, I_u] = BLS(I_l, I_u, flag(1+pad:row+pad,1+pad:col+pad), [row, col]);
    
    % Final Step for Output Image
    OutImg = RnS(I_l, I_u, flag(1+pad:row+pad,1+pad:col+pad), [row, col]);
   
end

function [l_pix, u_pix] = ABIC(pad, nImg, flag)
    win = 3;
    mid = pad+1;
    
    % Default to corrupted value in case the methods fails
    l_pix = nImg(mid,mid);
    u_pix = nImg(mid,mid);
    
    while win <= 2*pad+1
        rng = (win-1)/2;
        win = win + 2;
        temp = reshape(nImg(mid-rng:mid+rng,mid-rng:mid+rng).*flag(mid-rng:mid+rng,mid-rng:mid+rng), 1,[]);
        temp(temp == 0) = [];
        if ~isempty(temp)
            l_pix = min(temp);
            u_pix = max(temp);
            break;
        end
    end
end

% Bound Limit Step to smooth out the lower and upper bound images
function [O_l, O_u] = BLS(I_l, I_u, flag, size)

    row = size(1);
    col = size(2);

    O_l = zeros(size);
    O_u = zeros(size);

    % Padding for Image Smoothing
    pad = 1;
    I_l = padarray(I_l, [pad, pad]);
    I_u = padarray(I_u, [pad, pad]);
    
    k = [
        0 1 0;
        1 1 1;
        0 1 0
        ];
    k = k./sum(k,'all');
    
    for i = 1 + pad : row + pad
        for j = 1 + pad : col + pad
            if flag(i-pad,j-pad)==0
                O_l(i-pad,j-pad) = sum(I_l(i-pad:i+pad,j-pad:j+pad).*k,'all');
                O_u(i-pad,j-pad) = sum(I_u(i-pad:i+pad,j-pad:j+pad).*k,'all');
            else
                O_l(i-pad,j-pad) = I_l(i,j);
                O_u(i-pad,j-pad) = I_u(i,j);
            end
        end
    end
end

function out = RnS(I_1, I_2, flag, size)    

    row = size(1);
    col = size(2);
    out = uint8(zeros(size));
    
    % Image Padding
    pad = 1;
    I_1 = padarray(I_1, [pad, pad]);
    I_2 = padarray(I_2, [pad, pad]);
    
    k = [
        0 1 0;
        1 1 1;
        0 1 0
        ];
    k = k./sum(k,'all');
    
    % Recombination Step
    oImg = (I_1.* 0.5 + I_2.* 0.5);
    
    % Smoothing Step
    for i = 1 + pad : row + pad
        for j = 1 + pad : col + pad
            if flag(i-pad,j-pad)==0
                out(i-pad,j-pad) = uint8(sum(oImg(i-pad:i+pad,j-pad:j+pad).*k,'all'));
            else
                out(i-pad,j-pad) = uint8(oImg(i,j));
            end
        end
    end
end