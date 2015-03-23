function compare(ORIGINAL,OUTPUTED)
origImg = imread(ORIGINAL);
distImg = imread(OUTPUTED);
[PSNR_val,MSE] = PSNR(origImg,distImg);
[TOTAL_COEF,CHANGED,PERCENTAGE] = changed(origImg,distImg);
disp('PIXEL DOMAIN');
disp(['PSNR = ' num2str(mean(PSNR_val))]);
disp(['MSE = ' num2str(mean(MSE))]);
disp(['TOTAL BITS = ' num2str(TOTAL_COEF)]);
disp(['CHANGED BITS = ' num2str(CHANGED)]);
disp(['CHANGE BIT RATIO = ' num2str(PERCENTAGE)]);


%[ORIGINAL_DCT_A,ORIGINAL_DCT_B,ORIGINAL_DCT_C] = DCTBlocks(ORIGINAL);
%[OUTPUTED_DCT_A,OUTPUTED_DCT_B,OUTPUTED_DCT_C] = DCTBlocks(OUTPUTED);
%[TOTAL_COEF,CHANGED,PERCENTAGE] = count_dct_changes(ORIGINAL_DCT_A,ORIGINAL_DCT_B,ORIGINAL_DCT_C,OUTPUTED_DCT_A,OUTPUTED_DCT_B,OUTPUTED_DCT_C);
%disp('FREQ DOMAIN');
%disp(['TOTAL BITS = ' num2str(TOTAL_COEF)]);
%disp(['CHANGED BITS = ' num2str(CHANGED)]);
%disp(['CHANGE RATIO = ' num2str(PERCENTAGE)]);