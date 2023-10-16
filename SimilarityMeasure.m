%------------------ Similarity function -----------------------------------
function [distance] = SimilarityMeasure(FirstImage,SecondImage)
    distance = sum(abs((FirstImage-SecondImage)./(FirstImage+SecondImage+1e-10)));
end % end of Similarity function