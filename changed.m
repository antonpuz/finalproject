function [TOTAL_COEF,CHANGED,PERCENTAGE] = changed(origImg,distImg)
TOTAL_COEF = 0;
CHANGED = 0;
for d=1:3
    for j=1:512
        for i=1:512
            BITS = double(bitxor(origImg(i,j,d),distImg(i,j,d)));
            [nothing,e]=log2(max(BITS));
            CHANGED = CHANGED + sum(rem(floor(BITS*pow2(1-max(1,e):0)),2));
            TOTAL_COEF = TOTAL_COEF + 8;
        end
    end
end

PERCENTAGE = (CHANGED/TOTAL_COEF)*100;


