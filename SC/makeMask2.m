function [outIm] = makeMask2(bg, im, tol)
outIm = ((im > bg + tol) | (im < bg - tol));
end