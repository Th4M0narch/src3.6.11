function M = calculate_psd_metrics(expMass,qesMass,expNumber,qesNumber, ...
    d50Exp,d50Qes,sigmaExp,sigmaQes)
% Calculate PSD comparison metrics for one or more matched sampler rows.

expMass = as_matrix(expMass);
qesMass = as_matrix(qesMass);
expNumber = as_matrix(expNumber);
qesNumber = as_matrix(qesNumber);

n = max([size(expMass,1),size(qesMass,1), ...
         size(expNumber,1),size(qesNumber,1), ...
         numel(d50Exp),numel(d50Qes), ...
         numel(sigmaExp),numel(sigmaQes)]);

d50Exp = as_vector(d50Exp,n);
d50Qes = as_vector(d50Qes,n);
sigmaExp = as_vector(sigmaExp,n);
sigmaQes = as_vector(sigmaQes,n);

M = struct( ...
    'd50_error_relative',NaN(n,1), ...
    'sigma_error',NaN(n,1), ...
    'mass_RMSE',NaN(n,1), ...
    'number_RMSE',NaN(n,1), ...
    'JS_divergence',NaN(n,1));

for k = 1:n
    if k <= size(expMass,1) && k <= size(qesMass,1)
        pm = expMass(k,:);
        qm = qesMass(k,:);
        if numel(pm) == numel(qm) && any(isfinite(pm)) && any(isfinite(qm))
            M.mass_RMSE(k) = rmse_finite(pm,qm);
            M.JS_divergence(k) = jensen_shannon(pm,qm);
        end
    end
    if k <= size(expNumber,1) && k <= size(qesNumber,1)
        pn = expNumber(k,:);
        qn = qesNumber(k,:);
        if numel(pn) == numel(qn) && any(isfinite(pn)) && any(isfinite(qn))
            M.number_RMSE(k) = rmse_finite(pn,qn);
        end
    end
    if isfinite(d50Exp(k)) && d50Exp(k) > 0 && isfinite(d50Qes(k))
        M.d50_error_relative(k) = (d50Qes(k)-d50Exp(k)) / d50Exp(k);
    end
    if isfinite(sigmaExp(k)) && isfinite(sigmaQes(k))
        M.sigma_error(k) = sigmaQes(k) - sigmaExp(k);
    end
end
end

function x = as_matrix(x)
if isempty(x)
    x = zeros(0,0);
    return;
end
if isvector(x)
    x = x(:)';
end
end

function v = as_vector(v,n)
if isempty(v)
    v = NaN(n,1);
    return;
end
v = v(:);
if numel(v) < n
    v(end+1:n,1) = NaN;
end
end

function value = rmse_finite(a,b)
b = b(:)';
q = isfinite(a) & isfinite(b);
if ~any(q)
    value = NaN;
else
    value = sqrt(mean((a(q)-b(q)).^2));
end
end

function value = jensen_shannon(p,q)
p = p(:);
q = q(:);
qmask = isfinite(p) & isfinite(q) & p >= 0 & q >= 0;
p = p(qmask);
q = q(qmask);
if isempty(p) || sum(p) <= 0 || sum(q) <= 0
    value = NaN;
    return;
end
p = p / sum(p);
q = q / sum(q);
m = 0.5 * (p + q);
value = 0.5 * kl_divergence(p,m) + 0.5 * kl_divergence(q,m);
end

function value = kl_divergence(p,q)
qmask = p > 0 & q > 0;
if ~any(qmask)
    value = 0;
else
    value = sum(p(qmask) .* log(p(qmask)./q(qmask)));
end
end
