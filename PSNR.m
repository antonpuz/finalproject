         function [PSNR,MSE] = psnr(origImg,distImg)

            origImg = double(origImg);
            distImg = double(distImg);

            [M N D] = size(origImg);
            error = origImg - distImg;
            MSE = sum(sum(error .* error)) / (M * N * D);
            maxPixelValue = max(max(max(origImg)));

            PSNR = 10*log10(maxPixelValue*maxPixelValue/MSE);


        end