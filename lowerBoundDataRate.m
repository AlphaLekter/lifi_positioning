function dRate = lowerBoundDataRate(impResp, B, p, R_pd, q, N)
% Lower bound channel capacity
%   @ARTICLE{6636053,
%   author={Wang, Jun-Bo and Hu, Qing-Song and Wang, Jiangzhou and Chen, Ming and Wang, Jin-Yuan},
%   journal={Journal of Lightwave Technology},
%   title={Tight Bounds on Channel Capacity for Dimmable Visible Light Communications},
%   year={2013},
%   volume={31},
%   number={23},
%   pages={3771-3779},
%   doi={10.1109/JLT.2013.2286088}}
%
SNR = (exp(1)/2*pi)* ((p* R_pd * impResp)/(q*N*B));
dRate = B * log2(1+SNR);
end