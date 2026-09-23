function M = cesta_metrics(pred, obs)
% CESTA/QES model-performance metrics.
% pred = model prediction, obs = observation.
%
% FAC2 = fraction within factor 2
% FAC5 = fraction within factor 5
% FAC10 = fraction within factor 10
% MG   = exp(mean(log(obs/pred))) for positive pairs
% VG   = exp(mean((log(obs)-log(pred)).^2)) for positive pairs
% FB   = 2*(mean(obs)-mean(pred))/(mean(obs)+mean(pred))
% NMSE = mean((pred-obs).^2)/(mean(pred)*mean(obs))
% NAD  = mean(abs(pred-obs))/(mean(pred)+mean(obs))
% CC   = linear correlation coefficient

pred = pred(:); obs = obs(:);
q = isfinite(pred) & isfinite(obs) & pred >= 0 & obs >= 0;
pred = pred(q); obs = obs(q);

M = struct('N',0,'FAC2',NaN,'FAC5',NaN,'FAC10',NaN,'MG',NaN,'VG',NaN,...
           'FB',NaN,'NMSE',NaN,'NAD',NaN,'CC',NaN);
M.N = numel(pred);
if isempty(pred), return; end

pos = pred > 0 & obs > 0;
if any(pos)
    r = pred(pos)./obs(pos);
    M.FAC2 = mean(r >= 0.5 & r <= 2.0);
    M.FAC5 = mean(r >= 0.2 & r <= 5.0);
    M.FAC10 = mean(r >= 0.1 & r <= 10.0);
    dlog = log(obs(pos)) - log(pred(pos));
    M.MG = exp(mean(dlog));
    M.VG = exp(mean(dlog.^2));
end

mp = mean(pred); mo = mean(obs);
den = mp + mo;
if den ~= 0
    M.FB  = 2*(mo-mp)/den;
    M.NAD = mean(abs(pred-obs))/den;
end
if mp > 0 && mo > 0
    M.NMSE = mean((pred-obs).^2)/(mp*mo);
end
if numel(pred) >= 2 && std(pred) > 0 && std(obs) > 0
    C = corrcoef(pred,obs);
    M.CC = C(1,2);
end
end
