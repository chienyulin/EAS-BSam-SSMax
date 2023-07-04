function [ad_map, qdm_report] = EAS_BSam(HIM,GT,para_set)
% 
% EAS-BSam implementation
%
% INPUTS:
%   - HIM:      hyperspectral imagecube (rows by columns by bands)
%   - GT:       ground truth (rows by columns)
%   - para_set: parameters setting including maximum iteration number,
%               stopping rule criteria, and an initial random seed.
% OUTPUTS:
%   - ad_map:       final anomaly detection map (rows by columns)
%   - qdm_report:   quantitative detection measures report
%
% REFERENCE:
%   To-Be-Added
%
% Written and sorted by © 2023 Chien-Yu Lin.
%
fprintf('__initial__\n');
% random seed setting for reproducibility
rng(para_set.loop_rng_seed)
num_iter = 0;
HIM_feedback = [];
ad_map_stack = [];
b_map_stack = [];
t_stack = [];
pre_t_stack = [];
TI_stack = [];
AUCresult_stack = [];
currentmaps = [];
[m,n,~] = size(HIM);

FAP = 10^-4;
n_bs = NWHFC(HIM, FAP);
para_set.n_bs = n_bs;   % this number will be used in BSam as n_bs or in EAS as p

%% EAS space
HIM = rerange_data(HIM);
HIM = EAS(HIM,para_set);

%% initial state
% random select n_bs bands and reselect for each iteration
BSam_HIM = func_rand_band_selection(HIM,para_set);
% update random seed
para_set.loop_rng_seed = para_set.loop_rng_seed+1;

HIM_pro = BSam_HIM;
% normalize each band to -1 and 1
HIM_pro = rerange_data(HIM_pro);
% get AD map
t1start = tic;
% replace the line with your anomaly detector, for example, RXAD, PTA,
% RGAE,... etc, and produce the ad_map
ad_map = reshape(mat2gray(RxDetector(ToVector(HIM_pro)')),m,n);
t1end = toc(t1start);
% save current ad_map to ad_map_stack
ad_map_stack = cat(3,ad_map_stack,ad_map);

% calculate AUC for AD method without iteration
AUCresult = get_AUC_result(ad_map,GT);
AUCresult_stack = cat(1,AUCresult_stack,AUCresult);

% Ostu's method, save threhold from each iteration
t = graythresh(ad_map);
t_stack = cat(1,t_stack,t);

% binary map created and save every binary map from each iterations
% anomaly=1;background=0
b_map = imbinarize(ad_map,t);
b_map_stack = cat(3,b_map_stack,b_map);

% Spatial filter map created
sf_map = EPF(3,1,ad_map,double(b_map));

% max_map
max_map = max(ad_map,sf_map);

% save current maps to the buffer "currentmaps"
currentmaps.ad_map = ad_map;
currentmaps.sf_map = sf_map;
currentmaps.max_map = max_map;

% create feedback cube
HIM_feedback = cat(3,HIM_feedback,sf_map);

% add back to original band sampling cube
BSam_HIM = func_rand_band_selection(HIM,para_set);
% update random seed
para_set.loop_rng_seed = para_set.loop_rng_seed+1;
HIM_pro = cat(3,BSam_HIM,HIM_feedback);

% save first AD map to pre map
pre_map = currentmaps.ad_map;

clearvars currentmaps
%% iteration
while num_iter < para_set.max_iter
    num_iter = num_iter+1;
    fprintf('iteration No. %d', num_iter);

    % normalize each band to -1 and 1
    HIM_pro = rerange_data(HIM_pro);

    % anomaly detection methods, ad_map created
    ad_map = reshape(mat2gray(RxDetector(ToVector(HIM_pro)')),m,n);
    ad_map_stack = cat(3,ad_map_stack,ad_map);

    % AUC
    AUCresult = get_AUC_result(ad_map,GT);
    AUCresult_stack = cat(1,AUCresult_stack,AUCresult);

    %% examine stopping rule, TI
    % otsu's method
    [TI, b_map, ~, t, pre_t] = cal_ti_otsu(ad_map, pre_map);
    b_map_stack = cat(3,b_map_stack,b_map);
    pre_t_stack = cat(1,pre_t_stack,pre_t);
    TI_stack = cat(1,TI_stack,TI);
    t_stack = cat(1,t_stack,t);
    fprintf(', TI is %f \n', TI);

    % similarity greater than e, break the while loop
    % if similarity not greater than e, continue the while loop
    if TI > para_set.e
        break;
    end

    % apply the spatial filter on ad_map
    sf_map = EPF(3,1,ad_map,double(b_map));

    % max_map
    max_map = max(ad_map,sf_map);

    % save current maps to the buffer "currentmaps"
    currentmaps.ad_map = ad_map;
    currentmaps.sf_map = sf_map;
    currentmaps.max_map = max_map;

    % create feedback cube
    HIM_feedback = cat(3,HIM_feedback,sf_map);

    % prepare HIM for next iteration,
    BSam_HIM = func_rand_band_selection(HIM,para_set);
    % update random seed
    para_set.loop_rng_seed = para_set.loop_rng_seed+1;
    HIM_pro = cat(3,BSam_HIM,HIM_feedback);

    % save current AD map to pre map
    pre_map = currentmaps.ad_map;
end
% iteration computation time
t2end = toc(t1start);
fprintf('iteration finished.\n');

%% create Quantitative Detection Measures (QDM) report
j = (0:num_iter)';
TI = [NaN; TI_stack];
CPUTime = [t1end;nan(numel(j)-2,1);t2end];
report_matrix = [j,TI,AUCresult_stack,CPUTime];
report_title = ["iteration_no.","TI","df","adp","bdp","jad","jbs","adbs","snpr","oadp","cputime"];
qdm_report = array2table(report_matrix);
qdm_report.Properties.VariableNames = report_title;

end
