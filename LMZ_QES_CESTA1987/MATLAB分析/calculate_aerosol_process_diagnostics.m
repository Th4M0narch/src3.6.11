function D = calculate_aerosol_process_diagnostics(inputs)
% Calculate deposition, coagulation, mass-balance, and distance diagnostics.

D = struct();
D.particle = table();
D.distance = table();
D.coagulation = table();
D.massBalance = table();
D.massBalanceSummary = table();
D.depositionBins = table();
D.psd = build_psd_table(inputs.psd);

if ~inputs.hasMass
    warning('Aerosol process diagnostics require aerosol_mass_diagnostic.csv.');
    return;
end

M = build_mass_step_table(inputs.mass);
if isempty(M)
    warning('Aerosol mass diagnostic contains no usable rows.');
    return;
end
M = attach_psd_fields(M,D.psd);

D.particle = summarize_particles(M);
D.distance = summarize_distance(D.particle);
D.coagulation = summarize_coagulation(D.particle);
D.massBalance = summarize_mass_balance_by_time(M);
D.massBalanceSummary = summarize_mass_balance_overall(M,D.particle);
D.depositionBins = summarize_deposition_bins(inputs.deposition);
end

function P = build_psd_table(T)
P = table();
if isempty(T)
    return;
end

time = numeric_column(getcol(T,{'time','time_s'},false));
particleID = numeric_column(getcol(T,{'particleID','particle_id','id'},false));
distance = numeric_column(getcol(T,{'distance','distance_m'},false));
if isempty(time) || isempty(particleID)
    warning('aerosol_psd.csv is missing time or particleID.');
    return;
end
if isempty(distance)
    x = numeric_column(getcol(T,{'x','x_m'},false));
    y = numeric_column(getcol(T,{'y','y_m'},false));
    if ~isempty(x) && ~isempty(y)
        distance = sqrt(x.^2 + y.^2);
    else
        distance = NaN(size(time));
    end
end

mass = zeros(numel(time),9);
for b = 1:9
    value = numeric_column(getcol(T, ...
        {sprintf('bin%d_mass_fraction',b)},false));
    if ~isempty(value)
        mass(:,b) = value;
    end
end
d50 = numeric_column(getcol(T,{'d50_um','d50'},false));
sigma = numeric_column(getcol(T,{'sigma_g','geometric_sigma'},false));
if isempty(d50)
    d50 = NaN(size(time));
end
if isempty(sigma)
    sigma = NaN(size(time));
end

binDiameter = representative_diameter_um();
for k = 1:numel(time)
    if ~isfinite(d50(k))
        d50(k) = mass_median_diameter(binDiameter,mass(k,:));
    end
    if ~isfinite(sigma(k))
        sigma(k) = mass_geometric_sigma(binDiameter,mass(k,:));
    end
end

P = table(time,particleID,distance,d50,sigma, ...
    'VariableNames',{'time','particleID','distance','d50','sigma'});
end

function M = build_mass_step_table(T)
M = table();
time = numeric_column(getcol(T,{'time','time_s'},false));
particleID = numeric_column(getcol(T,{'particleID','particle_id','id'},false));
if isempty(time) || isempty(particleID)
    return;
end

n = numel(time);
before = optional_column(T,{'uo2f2_before_kg'},n,NaN);
afterCoagulation = optional_column(T,{'uo2f2_after_coagulation_kg'},n,NaN);
afterChemistry = optional_column(T,{'uo2f2_after_chemistry_kg'},n,NaN);
afterDeposition = optional_column(T,{'uo2f2_after_deposition_kg'},n,NaN);
deposited = optional_column(T,{'uo2f2_deposited_kg'},n,NaN);
stepError = optional_column(T,{'uo2f2_step_error_kg'},n,NaN);

coagMassBefore = optional_column(T,{'coagulation_mass_before_kg'},n,NaN);
coagMassAfter = optional_column(T,{'coagulation_mass_after_kg'},n,NaN);
coagMassError = optional_column(T,{'coagulation_mass_error_kg'},n,NaN);
coagNumberBefore = optional_column(T,{'coagulation_number_before'},n,NaN);
coagNumberAfter = optional_column(T,{'coagulation_number_after'},n,NaN);

fError = optional_column(T,{'F_error_kg'},n,NaN);
fTotal = optional_column(T,{'F_total_kg'},n,NaN);
waterBefore = optional_column(T,{'water_before_chemistry_kg'},n,NaN);
waterAfter = optional_column(T,{'water_after_deposition_kg'},n,NaN);
waterDeposited = optional_column(T,{'water_deposited_kg'},n,NaN);
waterError = optional_column(T,{'water_error_kg'},n,NaN);

stepError(isnan(stepError)) = afterDeposition(isnan(stepError)) ...
    + deposited(isnan(stepError)) - before(isnan(stepError));
coagMassError(isnan(coagMassError)) = ...
    coagMassAfter(isnan(coagMassError)) - coagMassBefore(isnan(coagMassError));

depositionFraction = safe_ratio(deposited,before);
coagRelativeError = safe_ratio(abs(coagMassError),coagMassBefore);
massStepRelativeError = safe_ratio(abs(stepError),before);

meanMassBefore = safe_ratio(coagMassBefore,coagNumberBefore);
meanMassAfter = safe_ratio(coagMassAfter,coagNumberAfter);
coagDiameterGrowth = (safe_ratio(meanMassAfter,meanMassBefore)).^(1.0/3.0);
coagDiameterGrowth(~isfinite(coagDiameterGrowth)) = NaN;

M = table(time,particleID,before,afterCoagulation,afterChemistry, ...
    afterDeposition,deposited,stepError, ...
    coagMassBefore,coagMassAfter,coagMassError,coagRelativeError, ...
    coagNumberBefore,coagNumberAfter,coagDiameterGrowth, ...
    depositionFraction,massStepRelativeError,fError,fTotal, ...
    waterBefore,waterAfter,waterDeposited,waterError, ...
    'VariableNames',{'time','particleID','uo2f2_before','uo2f2_after_coagulation', ...
    'uo2f2_after_chemistry','uo2f2_after_deposition','uo2f2_deposited', ...
    'uo2f2_step_error', ...
    'coagulation_mass_before','coagulation_mass_after', ...
    'coagulation_mass_error','coagulation_mass_relative_error', ...
    'coagulation_number_before','coagulation_number_after', ...
    'coagulation_diameter_growth','deposition_fraction', ...
    'uo2f2_step_relative_error','F_error','F_total', ...
    'water_before','water_after_deposition','water_deposited', ...
    'water_error'});

[~,order] = sortrows([M.particleID M.time],[1 2]);
M = M(order,:);
end

function M = attach_psd_fields(M,P)
n = height(M);
M.distance = NaN(n,1);
M.d50 = NaN(n,1);
M.sigma = NaN(n,1);
if isempty(P)
    return;
end

particleIDs = unique(M.particleID(isfinite(M.particleID)));
for k = 1:numel(particleIDs)
    pid = particleIDs(k);
    mRows = find(M.particleID == pid);
    pRows = find(P.particleID == pid);
    if isempty(pRows)
        continue;
    end

    pt = P.time(pRows);
    pd = P.distance(pRows);
    p50 = P.d50(pRows);
    ps = P.sigma(pRows);
    finiteTime = isfinite(pt);
    pt = pt(finiteTime);
    pd = pd(finiteTime);
    p50 = p50(finiteTime);
    ps = ps(finiteTime);
    if isempty(pt)
        continue;
    end
    [pt,uniqueIndex] = unique(pt,'last');
    pd = pd(uniqueIndex);
    p50 = p50(uniqueIndex);
    ps = ps(uniqueIndex);

    queryTime = M.time(mRows);
    if isscalar(pt)
        M.distance(mRows) = repmat(pd,1,numel(mRows));
        M.d50(mRows) = repmat(p50,1,numel(mRows));
        M.sigma(mRows) = repmat(ps,1,numel(mRows));
    else
        M.distance(mRows) = interp1(pt,pd,queryTime,'nearest','extrap');
        M.d50(mRows) = interp1(pt,p50,queryTime,'nearest','extrap');
        M.sigma(mRows) = interp1(pt,ps,queryTime,'nearest','extrap');
    end
end
end

function S = summarize_particles(M)
particleIDs = unique(M.particleID(isfinite(M.particleID)));
n = numel(particleIDs);
particleID = NaN(n,1);
time_first_s = NaN(n,1);
time_last_s = NaN(n,1);
distance_m = NaN(n,1);
d50_initial_um = NaN(n,1);
d50_final_um = NaN(n,1);
d50_growth = NaN(n,1);
coagulation_growth = NaN(n,1);
number_survival = NaN(n,1);
deposited_mass_kg = NaN(n,1);
airborne_mass_kg = NaN(n,1);
deposition_fraction = NaN(n,1);
mass_loss_fraction = NaN(n,1);
airborne_fraction = NaN(n,1);
steps = NaN(n,1);

for k = 1:n
    pid = particleIDs(k);
    rows = find(M.particleID == pid);
    particleID(k) = pid;
    steps(k) = numel(rows);
    times = M.time(rows);
    qTime = isfinite(times);
    if any(qTime)
        time_first_s(k) = min(times(qTime));
        time_last_s(k) = max(times(qTime));
    end

    lastDistance = M.distance(rows);
    qDistance = isfinite(lastDistance);
    if any(qDistance)
        distance_m(k) = lastDistance(find(qDistance,1,'last'));
    end

    d50Values = M.d50(rows);
    qD50 = isfinite(d50Values) & d50Values > 0;
    if any(qD50)
        d50_initial_um(k) = d50Values(find(qD50,1,'first'));
        d50_final_um(k) = d50Values(find(qD50,1,'last'));
        d50_growth(k) = d50_final_um(k) / d50_initial_um(k);
    end

    growth = M.coagulation_diameter_growth(rows);
    qGrowth = isfinite(growth) & growth > 0;
    if any(qGrowth)
        coagulation_growth(k) = prod(growth(qGrowth));
    else
        coagulation_growth(k) = 1.0;
    end

    ratios = safe_ratio(M.coagulation_number_after(rows), ...
        M.coagulation_number_before(rows));
    qRatio = isfinite(ratios) & ratios >= 0;
    if any(qRatio)
        number_survival(k) = prod(ratios(qRatio));
    else
        number_survival(k) = 1.0;
    end

    airborne = M.uo2f2_after_deposition(rows);
    qAirborne = isfinite(airborne);
    if any(qAirborne)
        airborne_mass_kg(k) = airborne(find(qAirborne,1,'last'));
    else
        airborne_mass_kg(k) = 0.0;
    end
    deposited_mass_kg(k) = sum_finite(M.uo2f2_deposited(rows));
    total = airborne_mass_kg(k) + deposited_mass_kg(k);
    if total > 0
        deposited_mass_kg(k) = max(0.0,deposited_mass_kg(k));
        deposition_fraction(k) = deposited_mass_kg(k) / total;
        mass_loss_fraction(k) = deposition_fraction(k);
        airborne_fraction(k) = airborne_mass_kg(k) / total;
    end
end

S = table(particleID,time_first_s,time_last_s,distance_m, ...
    d50_initial_um,d50_final_um,d50_growth,coagulation_growth, ...
    number_survival,deposited_mass_kg,airborne_mass_kg, ...
    deposition_fraction,mass_loss_fraction,airborne_fraction,steps);
end

function T = summarize_distance(S)
T = table();
if isempty(S)
    return;
end
edges = [0 10 20 40 70 100 200 500 Inf];
distance_lower_m = NaN(numel(edges)-1,1);
distance_upper_m = NaN(numel(edges)-1,1);
samples = NaN(numel(edges)-1,1);
d50_initial_um = NaN(numel(edges)-1,1);
d50_final_um = NaN(numel(edges)-1,1);
d50_growth = NaN(numel(edges)-1,1);
coagulation_growth = NaN(numel(edges)-1,1);
deposition_fraction = NaN(numel(edges)-1,1);
mass_loss_fraction = NaN(numel(edges)-1,1);
deposited_mass_kg = NaN(numel(edges)-1,1);
airborne_mass_kg = NaN(numel(edges)-1,1);

d = S.distance_m;
for k = 1:numel(edges)-1
    if isinf(edges(k+1))
        q = isfinite(d) & d >= edges(k);
    else
        q = isfinite(d) & d >= edges(k) & d < edges(k+1);
    end
    distance_lower_m(k) = edges(k);
    distance_upper_m(k) = edges(k+1);
    samples(k) = sum(q);
    if ~any(q)
        continue;
    end
    d50_initial_um(k) = mean_finite(S.d50_initial_um(q));
    d50_final_um(k) = mean_finite(S.d50_final_um(q));
    d50_growth(k) = mean_finite(S.d50_growth(q));
    coagulation_growth(k) = mean_finite(S.coagulation_growth(q));
    deposition_fraction(k) = mean_finite(S.deposition_fraction(q));
    mass_loss_fraction(k) = mean_finite(S.mass_loss_fraction(q));
    deposited_mass_kg(k) = sum_finite(S.deposited_mass_kg(q));
    airborne_mass_kg(k) = sum_finite(S.airborne_mass_kg(q));
end

T = table(distance_lower_m,distance_upper_m,samples, ...
    d50_initial_um,d50_final_um,d50_growth,coagulation_growth, ...
    deposition_fraction,mass_loss_fraction,deposited_mass_kg, ...
    airborne_mass_kg);
end

function T = summarize_coagulation(S)
T = table();
if isempty(S)
    return;
end
tracked_particles = height(S);
mean_d50_growth = mean_finite(S.d50_growth);
median_d50_growth = median_finite(S.d50_growth);
mean_coagulation_growth = mean_finite(S.coagulation_growth);
median_coagulation_growth = median_finite(S.coagulation_growth);
mean_number_survival = mean_finite(S.number_survival);
median_number_survival = median_finite(S.number_survival);
mean_diameter_increase_percent = 100.0 * (mean_coagulation_growth - 1.0);
T = table(tracked_particles,mean_d50_growth,median_d50_growth, ...
    mean_coagulation_growth,median_coagulation_growth, ...
    mean_number_survival,median_number_survival, ...
    mean_diameter_increase_percent);
end

function T = summarize_mass_balance_by_time(M)
times = unique(M.time(isfinite(M.time)));
n = numel(times);
time_s = NaN(n,1);
particles = NaN(n,1);
uo2f2_before_kg = NaN(n,1);
uo2f2_after_coagulation_kg = NaN(n,1);
uo2f2_after_chemistry_kg = NaN(n,1);
uo2f2_after_deposition_kg = NaN(n,1);
uo2f2_deposited_kg = NaN(n,1);
uo2f2_step_error_kg = NaN(n,1);
uo2f2_step_relative_error = NaN(n,1);
coagulation_mass_error_kg = NaN(n,1);
coagulation_relative_error = NaN(n,1);
number_before = NaN(n,1);
number_after = NaN(n,1);
number_ratio = NaN(n,1);
F_error_kg = NaN(n,1);
F_relative_error = NaN(n,1);
water_error_kg = NaN(n,1);
water_relative_error = NaN(n,1);
closure_pass = false(n,1);

for k = 1:n
    q = M.time == times(k);
    time_s(k) = times(k);
    particles(k) = sum(q);
    uo2f2_before_kg(k) = sum_finite(M.uo2f2_before(q));
    uo2f2_after_coagulation_kg(k) = ...
        sum_finite(M.uo2f2_after_coagulation(q));
    uo2f2_after_chemistry_kg(k) = ...
        sum_finite(M.uo2f2_after_chemistry(q));
    uo2f2_after_deposition_kg(k) = ...
        sum_finite(M.uo2f2_after_deposition(q));
    uo2f2_deposited_kg(k) = sum_finite(M.uo2f2_deposited(q));
    uo2f2_step_error_kg(k) = sum_finite(M.uo2f2_step_error(q));
    uo2f2_step_relative_error(k) = abs_relative( ...
        uo2f2_step_error_kg(k),uo2f2_before_kg(k));
    coagulation_mass_error_kg(k) = ...
        sum_finite(M.coagulation_mass_error(q));
    coagulation_relative_error(k) = abs_relative( ...
        coagulation_mass_error_kg(k),sum_finite(M.coagulation_mass_before(q)));
    number_before(k) = sum_finite(M.coagulation_number_before(q));
    number_after(k) = sum_finite(M.coagulation_number_after(q));
    number_ratio(k) = safe_ratio(number_after(k),number_before(k));
    F_error_kg(k) = sum_finite(M.F_error(q));
    F_expected = sum_finite(M.F_total(q)) - F_error_kg(k);
    F_relative_error(k) = abs_relative(F_error_kg(k),F_expected);
    water_error_kg(k) = sum_finite(M.water_error(q));
    water_relative_error(k) = abs_relative( ...
        water_error_kg(k),sum_finite(M.water_before(q)));
    closure_pass(k) = closure_ok(uo2f2_step_error_kg(k),uo2f2_before_kg(k)) ...
        && closure_ok(coagulation_mass_error_kg(k), ...
            sum_finite(M.coagulation_mass_before(q))) ...
        && closure_ok(F_error_kg(k),F_expected) ...
        && closure_ok(water_error_kg(k),sum_finite(M.water_before(q)));
end

T = table(time_s,particles,uo2f2_before_kg, ...
    uo2f2_after_coagulation_kg,uo2f2_after_chemistry_kg, ...
    uo2f2_after_deposition_kg,uo2f2_deposited_kg, ...
    uo2f2_step_error_kg,uo2f2_step_relative_error, ...
    coagulation_mass_error_kg,coagulation_relative_error, ...
    number_before,number_after,number_ratio,F_error_kg,F_relative_error, ...
    water_error_kg,water_relative_error,closure_pass);
end

function T = summarize_mass_balance_overall(M,S)
uo2f2_deposited_kg = sum_finite(M.uo2f2_deposited);
uo2f2_airborne_final_kg = sum_finite(S.airborne_mass_kg);
uo2f2_total_kg = uo2f2_deposited_kg + uo2f2_airborne_final_kg;
uo2f2_deposition_fraction = safe_ratio( ...
    uo2f2_deposited_kg,uo2f2_total_kg);
uo2f2_step_error_kg = sum_finite(M.uo2f2_step_error);
uo2f2_step_relative_error = abs_relative( ...
    uo2f2_step_error_kg,sum_finite(M.uo2f2_before));
coagulation_mass_error_kg = sum_finite(M.coagulation_mass_error);
coagulation_relative_error = abs_relative( ...
    coagulation_mass_error_kg,sum_finite(M.coagulation_mass_before));
F_error_kg = sum_finite(M.F_error);
F_expected_kg = sum_finite(M.F_total) - F_error_kg;
F_relative_error = abs_relative(F_error_kg,F_expected_kg);
water_error_kg = sum_finite(M.water_error);
water_relative_error = abs_relative( ...
    water_error_kg,sum_finite(M.water_before));
T = table(height(S),height(M),uo2f2_airborne_final_kg, ...
    uo2f2_deposited_kg,uo2f2_total_kg,uo2f2_deposition_fraction, ...
    uo2f2_step_error_kg,uo2f2_step_relative_error, ...
    coagulation_mass_error_kg,coagulation_relative_error, ...
    F_error_kg,F_expected_kg,F_relative_error, ...
    water_error_kg,water_relative_error, ...
    'VariableNames',{'tracked_particles','diagnostic_rows', ...
    'uo2f2_airborne_final_kg','uo2f2_deposited_kg','uo2f2_total_kg', ...
    'uo2f2_deposition_fraction','uo2f2_step_error_kg', ...
    'uo2f2_step_relative_error','coagulation_mass_error_kg', ...
    'coagulation_relative_error','F_error_kg','F_expected_kg', ...
    'F_relative_error','water_error_kg','water_relative_error'});
end

function T = summarize_deposition_bins(Tdep)
T = table();
if isempty(Tdep)
    return;
end
time = numeric_column(getcol(Tdep,{'time','time_s'},false));
bin = numeric_column(getcol(Tdep,{'bin','bin_id'},false));
uo2f2 = numeric_column(getcol(Tdep,{'uo2f2_deposited_kg'},false));
water = numeric_column(getcol(Tdep,{'bound_water_deposited_kg'},false));
diameter = numeric_column(getcol(Tdep,{'diameter_um'},false));
reynolds = numeric_column(getcol(Tdep,{'reynolds_number'},false));
velocity = numeric_column(getcol(Tdep,{'settling_velocity_m_s'},false));
drag = numeric_column(getcol(Tdep,{'drag_coefficient'},false));
if isempty(bin) || isempty(uo2f2)
    warning('aerosol_deposition.csv is missing bin or deposited mass columns.');
    return;
end
if any(bin == 0)
    bin = bin + 1;
end
if isempty(water), water = NaN(size(uo2f2)); end
if isempty(diameter), diameter = NaN(size(uo2f2)); end
if isempty(reynolds), reynolds = NaN(size(uo2f2)); end
if isempty(velocity), velocity = NaN(size(uo2f2)); end
if isempty(drag), drag = NaN(size(uo2f2)); end
if isempty(time), time = NaN(size(uo2f2)); end

total = sum_finite(uo2f2);
bins = unique(bin(isfinite(bin)));
n = numel(bins);
bin_id = NaN(n,1);
events = NaN(n,1);
uo2f2_deposited_kg = NaN(n,1);
bound_water_deposited_kg = NaN(n,1);
deposition_mass_fraction = NaN(n,1);
mean_diameter_um = NaN(n,1);
mass_weighted_reynolds = NaN(n,1);
mass_weighted_settling_velocity_m_s = NaN(n,1);
mass_weighted_drag_coefficient = NaN(n,1);
first_time_s = NaN(n,1);
last_time_s = NaN(n,1);

for k = 1:n
    q = bin == bins(k);
    bin_id(k) = bins(k);
    events(k) = sum(q);
    uo2f2_deposited_kg(k) = sum_finite(uo2f2(q));
    bound_water_deposited_kg(k) = sum_finite(water(q));
    deposition_mass_fraction(k) = safe_ratio( ...
        uo2f2_deposited_kg(k),total);
    mean_diameter_um(k) = mean_finite(diameter(q));
    mass_weighted_reynolds(k) = weighted_mean(reynolds(q),uo2f2(q));
    mass_weighted_settling_velocity_m_s(k) = ...
        weighted_mean(velocity(q),uo2f2(q));
    mass_weighted_drag_coefficient(k) = ...
        weighted_mean(drag(q),uo2f2(q));
    qt = q & isfinite(time);
    if any(qt)
        first_time_s(k) = min(time(qt));
        last_time_s(k) = max(time(qt));
    end
end

T = table(bin_id,events,uo2f2_deposited_kg, ...
    bound_water_deposited_kg,deposition_mass_fraction, ...
    mean_diameter_um,mass_weighted_reynolds, ...
    mass_weighted_settling_velocity_m_s, ...
    mass_weighted_drag_coefficient,first_time_s,last_time_s);
end

function v = optional_column(T,names,n,defaultValue)
v = numeric_column(getcol(T,names,false));
if isempty(v)
    v = repmat(defaultValue,n,1);
else
    v = v(:);
    if numel(v) < n
        v(end+1:n,1) = defaultValue;
    end
end
end

function v = numeric_column(v)
if isempty(v)
    return;
end
if ~isnumeric(v)
    v = str2double(to_cellstr_compat(v));
end
v = double(v(:));
end

function value = safe_ratio(numerator,denominator)
value = NaN(size(numerator));
q = isfinite(numerator) & isfinite(denominator) & denominator ~= 0;
value(q) = numerator(q) ./ denominator(q);
end

function value = sum_finite(v)
v = v(isfinite(v));
if isempty(v)
    value = 0.0;
else
    value = sum(v);
end
end

function value = mean_finite(v)
v = v(isfinite(v));
if isempty(v)
    value = NaN;
else
    value = mean(v);
end
end

function value = median_finite(v)
v = v(isfinite(v));
if isempty(v)
    value = NaN;
else
    value = median(v);
end
end

function value = weighted_mean(v,w)
q = isfinite(v) & isfinite(w) & w >= 0;
if ~any(q) || sum(w(q)) <= 0
    value = NaN;
else
    value = sum(v(q).*w(q)) / sum(w(q));
end
end

function value = abs_relative(errorValue,scaleValue)
if ~isfinite(errorValue) || ~isfinite(scaleValue)
    value = NaN;
else
    value = abs(errorValue) / max(abs(scaleValue),1.0e-30);
end
end

function pass = closure_ok(errorValue,scaleValue)
if ~isfinite(errorValue)
    pass = false;
else
    tolerance = max(1.0e-15,1.0e-10*max(abs(scaleValue),0.0));
    pass = abs(errorValue) <= tolerance;
end
end

function d = representative_diameter_um()
d = [0.2 sqrt(0.4*0.7) sqrt(0.7*1.1) sqrt(1.1*2.1) ...
    sqrt(2.1*3.3) sqrt(3.3*4.7) sqrt(4.7*5.8) sqrt(5.8*9.0) 12.0];
end

function d50 = mass_median_diameter(diameter_um,mass_fraction)
d50 = NaN;
if isempty(diameter_um) || sum(mass_fraction) <= 0
    return;
end
p = mass_fraction(:) / sum(mass_fraction);
c = cumsum(p);
idx = find(c >= 0.5,1,'first');
if isempty(idx)
    return;
end
if idx == 1
    d50 = diameter_um(1);
    return;
end
c0 = c(idx-1);
if p(idx) <= 0
    d50 = diameter_um(idx);
else
    d50 = diameter_um(idx-1) ...
        + (0.5-c0)/p(idx) * (diameter_um(idx)-diameter_um(idx-1));
end
end

function sigma = mass_geometric_sigma(diameter_um,mass_fraction)
sigma = NaN;
q = isfinite(diameter_um) & isfinite(mass_fraction) ...
    & diameter_um > 0 & mass_fraction > 0;
if ~any(q)
    return;
end
p = mass_fraction(q) / sum(mass_fraction(q));
meanLog = sum(p .* log(diameter_um(q)));
sigma = exp(sqrt(sum(p .* (log(diameter_um(q))-meanLog).^2)));
end
