function M = calculate_aerosol_psd_metrics(PSD,reference)
% Calculate PSD evolution and shape metrics against an initial/reference PSD.

n = numel(PSD.qes);
M = struct( ...
    'distance_m',NaN(n,1), ...
    'height_m',NaN(n,1), ...
    'time_s',NaN(n,1), ...
    'd50_um',NaN(n,1), ...
    'sigma_g',NaN(n,1), ...
    'd50_change',NaN(n,1), ...
    'sigma_change',NaN(n,1), ...
    'mass_RMSE',NaN(n,1), ...
    'number_RMSE',NaN(n,1), ...
    'JS_divergence',NaN(n,1));

if n == 0
    return;
end

initialD50 = qes_initial_value(PSD.qes,'d50');
initialSigma = qes_initial_value(PSD.qes,'sigma_g');
referenceMass = get_reference(reference, ...
    'mass_fraction_9bin_by_sampler','mass_fraction_by_sampler');
referenceNumber = get_reference(reference, ...
    'number_fraction_9bin_by_sampler','number_fraction_by_sampler');
hasReference = ~isempty(referenceMass) && ~isempty(referenceNumber);

for i = 1:n
    M.distance_m(i) = PSD.qes(i).distance;
    M.height_m(i) = PSD.qes(i).height;
    M.time_s(i) = PSD.qes(i).time;
    M.d50_um(i) = PSD.qes(i).d50;
    M.sigma_g(i) = PSD.qes(i).sigma_g;

    if isfinite(initialD50) && initialD50 > 0 && isfinite(M.d50_um(i))
        M.d50_change(i) = (M.d50_um(i)-initialD50) / initialD50;
    end
    if isfinite(initialSigma) && isfinite(M.sigma_g(i))
        M.sigma_change(i) = M.sigma_g(i) - initialSigma;
    end

    if hasReference
        p = referenceMass(:)';
        q = PSD.qes(i).mass_fraction(:)';
        if numel(p) == numel(q)
            M.mass_RMSE(i) = rmse_finite(p,q);
            M.JS_divergence(i) = jensen_shannon(p,q);
        end

        pn = referenceNumber(:)';
        qn = PSD.qes(i).number_fraction(:)';
        if numel(pn) == numel(qn)
            M.number_RMSE(i) = rmse_finite(pn,qn);
        end
    end
end
end

function value = qes_initial_value(qes,fieldName)
value = NaN;
if isempty(qes) || ~isfield(qes,fieldName)
    return;
end
times = [qes.time];
values = [qes.(fieldName)];
qPositive = isfinite(times) & isfinite(values) & values > 0;
qFinite = isfinite(times) & isfinite(values);
if any(qPositive)
    initialTime = min(times(qPositive));
    selected = qPositive ...
        & abs(times-initialTime) <= max(1.0e-12,1.0e-9*abs(initialTime));
    value = mean(values(selected));
elseif any(qFinite)
    initialTime = min(times(qFinite));
    selected = qFinite ...
        & abs(times-initialTime) <= max(1.0e-12,1.0e-9*abs(initialTime));
    value = mean(values(selected));
end
end

function values = get_reference(reference,preferredName,fallbackName)
values = [];
if ~isfield(reference,'exp')
    return;
end
if isfield(reference.exp,preferredName) ...
        && ~isempty(reference.exp.(preferredName)) ...
        && any(reference.exp.(preferredName) > 0)
    values = reference.exp.(preferredName);
elseif isfield(reference.exp,fallbackName) ...
        && numel(reference.exp.(fallbackName)) == 9 ...
        && any(reference.exp.(fallbackName) > 0)
    values = reference.exp.(fallbackName);
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
