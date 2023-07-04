function  [TI, b_current_map, b_pre_map, current_t, pre_t] = cal_ti_otsu(current_ad_map, pre_ad_map)
% 
% calculate stopping rule
%
% INPUTS:
%   - current_ad_map:   current anomaly detection map
%   - pre_ad_map:       anomaly detection map in the previous loop
% OUTPUTS:
%   - TI:               Tanimoto index
%   - b_current_map:    binary map of current_ad_map
%   - b_pre_map:        binary map of pre_ad_map
%   - current_t:        threshold of current_ad_map
%   - pre_t:            threshold of pre_ad_map
%
% Written and sorted by © 2023 Chien-Yu Lin.

current_t = graythresh(current_ad_map);
b_current_map = imbinarize(current_ad_map,current_t);

pre_t = graythresh(pre_ad_map);
b_pre_map = imbinarize(pre_ad_map,pre_t);

s1 = find(ToVector(b_current_map) == 1);
s2 = find(ToVector(b_pre_map) == 1);

TI = length(intersect(s1,s2))/ length(union(s1,s2));
