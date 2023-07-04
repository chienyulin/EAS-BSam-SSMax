function AUCresult = get_AUC_result(ad_map,GT)
% 
% create quantitative detection measures
%
% INPUTS:
%   - ad_map:   anomaly detection map (rows by columns)
%   - GT:       ground truth (rows by columns)
% OUTPUTS:
%   - AUCresult:    quantitative detection measures
%
% Written and sorted by © 2023 Chien-Yu Lin.
%
[~,AUCnor] = cal_AUC(ToVector(ad_map),GT,1,1);
% (D,F),ADP ((D,tau)), BDP (1-(F,tau), JTD, JBS, ADBS, SNPR, OADP

eva_cell.df     = AUCnor.PFPD;  % DF
eva_cell.adp    = AUCnor.tauPD;  % ADP = Dtau
eva_cell.jad    = AUCnor.PFPD + AUCnor.tauPD;  % JAD = DF+ADP = DF+Dtau
eva_cell.bdp    = 1 - AUCnor.tauPF;  % BDP = 1-Ftau
eva_cell.jbs    = AUCnor.PFPD + 1 - AUCnor.tauPF;  % jbs = DF+BDP = DF+(1-Ftau)
eva_cell.adbs   = AUCnor.tauPD - AUCnor.tauPF;  % ADBS = Dtau-Ftau
eva_cell.oadp   = AUCnor.tauPD + 1 - AUCnor.tauPF;  % OADP = ADP+BDP = Dtau+(1-Ftau)
eva_cell.snpr   = AUCnor.tauPD/AUCnor.tauPF;  % SNPR = Dtau/Ftau

AUCresult = round([eva_cell.df,eva_cell.adp,eva_cell.bdp,eva_cell.jad,eva_cell.jbs,eva_cell.adbs,eva_cell.snpr,eva_cell.oadp],4);

