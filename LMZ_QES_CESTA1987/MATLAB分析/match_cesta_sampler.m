function [matched, report] = match_cesta_sampler(expSamplers, qesSamplers, ...
    azimuthWeight, timeWeight)
% Match experimental CESTA samplers to QES sampler rows.
%
% Priority:
%   1) exact sampler_id
%   2) time window
%   3) azimuth-weighted nearest distance/height/time
%
% Inputs may be tables or struct arrays.  Missing azimuth metadata falls
% back to distance/height/time matching.  One QES row is used at most once.

if nargin < 3 || isempty(azimuthWeight)
    azimuthWeight = 10.0;
end
if nargin < 4 || isempty(timeWeight)
    timeWeight = 1.0e-3;
end
if ~isfinite(azimuthWeight) || azimuthWeight < 1.0
    error('azimuthWeight must be finite and >= 1.');
end
if ~isfinite(timeWeight) || timeWeight < 0.0
    error('timeWeight must be finite and nonnegative.');
end

E = normalize_samplers(expSamplers);
Q = normalize_samplers(qesSamplers);
nE = numel(E.id);
nQ = numel(Q.id);

qesRow = NaN(nE,1);
quality = cell(nE,1);
distanceDifference = NaN(nE,1);
heightDifference = NaN(nE,1);
azimuthDifference = NaN(nE,1);
timeDifference = NaN(nE,1);
matchingCost = NaN(nE,1);
usedQes = false(nQ,1);

for i = 1:nE
    selected = 0;
    selectedCost = Inf;

    if ~isempty(E.id{i}) && any(~cellfun(@isempty,Q.id))
        idCandidates = find(~usedQes & strcmp(E.id{i},Q.id));
        if ~isempty(idCandidates)
            [selected,selectedCost] = best_candidate(E,i,Q,idCandidates, ...
                azimuthWeight,timeWeight,true);
            if selected > 0
                quality{i} = 'sampler_id';
            end
        end
    end

    if selected == 0
        candidates = find(~usedQes);
        [selected,selectedCost,qualityText] = best_candidate(E,i,Q, ...
            candidates,azimuthWeight,timeWeight,true);
        quality{i} = qualityText;
    end

    if selected == 0
        quality{i} = 'none';
        continue;
    end

    usedQes(selected) = true;
    qesRow(i) = selected;
    [distanceDifference(i),heightDifference(i),azimuthDifference(i), ...
        timeDifference(i),~] = candidate_errors(E,i,Q,selected, ...
        azimuthWeight,timeWeight);
    matchingCost(i) = selectedCost;
end

experimentSamplerId = E.id;
qesSamplerId = cell(nE,1);
for i = 1:nE
    if isfinite(qesRow(i))
        qesSamplerId{i} = Q.id{qesRow(i)};
    else
        qesSamplerId{i} = '';
    end
end

report = table(experimentSamplerId,qesSamplerId, ...
    distanceDifference,heightDifference,azimuthDifference, ...
    timeDifference,quality,matchingCost, ...
    'VariableNames', ...
    {'experiment_sampler_id','qes_sampler_id', ...
     'distance_difference','height_difference','azimuth_difference', ...
     'time_difference','matching_quality','matching_cost'});

experimentRow = (1:nE)';
qesRowOut = qesRow;
experimentDistance = E.distance;
experimentHeight = E.height;
experimentAzimuth = E.azimuth;
experimentTime = E.time;
qesDistance = NaN(nE,1);
qesHeight = NaN(nE,1);
qesAzimuth = NaN(nE,1);
qesTime = NaN(nE,1);
for i = 1:nE
    if isfinite(qesRow(i))
        j = qesRow(i);
        qesDistance(i) = Q.distance(j);
        qesHeight(i) = Q.height(j);
        qesAzimuth(i) = Q.azimuth(j);
        qesTime(i) = Q.time(j);
    end
end

matched = table(experimentSamplerId,qesSamplerId,experimentRow,qesRowOut, ...
    experimentDistance,experimentHeight,experimentAzimuth,experimentTime, ...
    qesDistance,qesHeight,qesAzimuth,qesTime,quality,matchingCost, ...
    'VariableNames', ...
    {'experiment_sampler_id','qes_sampler_id','experiment_row','qes_row', ...
     'experiment_distance_m','experiment_height_m', ...
     'experiment_azimuth_deg','experiment_time_s', ...
     'qes_distance_m','qes_height_m','qes_azimuth_deg','qes_time_s', ...
     'matching_quality','matching_cost'});
end

function [selected,cost,quality] = best_candidate(E,index,Q,candidates, ...
    azimuthWeight,timeWeight,allowTimeWindowPriority)
selected = 0;
cost = Inf;
quality = 'position_time';
if isempty(candidates)
    quality = 'none';
    return;
end

active = candidates(:);
if allowTimeWindowPriority ...
        && isfinite(E.windowStart(index)) && isfinite(E.windowEnd(index))
    hasWindow = isfinite(Q.windowStart(active)) & isfinite(Q.windowEnd(active));
    if any(hasWindow)
        % A cumulative result at t > windowEnd is valid when its actual
        % sampling window matches. A time merely inside the window is not.
        sameWindow = hasWindow ...
            & abs(Q.windowStart(active)-E.windowStart(index)) <= 1e-8 ...
            & abs(Q.windowEnd(active)-E.windowEnd(index)) <= 1e-8;
        active = active(sameWindow);
        if isempty(active), return; end
        knownComplete = isfinite(Q.sampleComplete(active));
        if any(knownComplete)
            active = active(knownComplete & Q.sampleComplete(active) == 1);
            if isempty(active), return; end
        end
        quality = 'position_complete_window';
    else
        % Legacy non-PSD callers without window metadata retain time matching.
        time = Q.time(active);
        inWindow = isfinite(time) ...
            & time >= E.windowStart(index)-1e-12 ...
            & time <= E.windowEnd(index)+1e-12;
        if any(inWindow)
            active = active(inWindow);
            quality = 'position_time_window';
        end
    end
end

costs = NaN(numel(active),1);
for k = 1:numel(active)
    [~,~,~,~,costs(k)] = candidate_errors(E,index,Q,active(k), ...
        azimuthWeight,timeWeight);
end
[minimumCost,j] = min(costs);
if isfinite(minimumCost)
    % PSD sampler rows are cumulative in time.  If several rows are equally
    % good (for example, repeated output times for one sampler_id), use the
    % latest time within the selected candidate set.
    tied = find(abs(costs - minimumCost) ...
        <= max(1.0e-12,1.0e-12 * abs(minimumCost)));
    tiedTimes = Q.time(active(tied));
    finiteTime = isfinite(tiedTimes);
    if any(finiteTime)
        [~,latest] = max(tiedTimes(finiteTime));
        finiteRows = tied(finiteTime);
        selected = active(finiteRows(latest));
    else
        selected = active(tied(1));
    end
    cost = minimumCost;
end
end

function [distanceError,heightError,azimuthError,timeError,cost] = ...
    candidate_errors(E,index,Q,j,azimuthWeight,timeWeight)
distanceError = paired_absolute_error(E.distance(index),Q.distance(j));
heightError = paired_absolute_error(E.height(index),Q.height(j));
azimuthError = angular_difference(E.azimuth(index),Q.azimuth(j));
if isfinite(E.windowStart(index)) && isfinite(E.windowEnd(index)) ...
        && isfinite(Q.windowStart(j)) && isfinite(Q.windowEnd(j))
    timeError = abs(E.windowStart(index)-Q.windowStart(j)) ...
        + abs(E.windowEnd(index)-Q.windowEnd(j));
else
    timeError = time_window_error(E,index,Q.time(j));
end

if ~isfinite(distanceError)
    distanceError = 0.0;
end
if ~isfinite(heightError)
    heightError = 0.0;
end
if ~isfinite(azimuthError)
    azimuthError = 0.0;
end
if ~isfinite(timeError)
    timeError = 0.0;
end

hasEvidence = isfinite(E.distance(index)) && isfinite(Q.distance(j)) ...
    || isfinite(E.height(index)) && isfinite(Q.height(j)) ...
    || isfinite(E.azimuth(index)) && isfinite(Q.azimuth(j)) ...
    || isfinite(Q.time(j));
if hasEvidence
    cost = distanceError + heightError + azimuthWeight * azimuthError ...
        + timeWeight * timeError;
else
    cost = Inf;
end
end

function value = paired_absolute_error(first,second)
if isfinite(first) && isfinite(second)
    value = abs(first - second);
else
    value = NaN;
end
end

function value = angular_difference(first,second)
if isfinite(first) && isfinite(second)
    value = abs(mod(first - second + 180.0,360.0) - 180.0);
else
    value = NaN;
end
end

function value = time_window_error(E,index,qesTime)
if ~isfinite(qesTime)
    value = NaN;
elseif isfinite(E.windowStart(index)) && isfinite(E.windowEnd(index))
    if qesTime < E.windowStart(index)
        value = E.windowStart(index) - qesTime;
    elseif qesTime > E.windowEnd(index)
        value = qesTime - E.windowEnd(index);
    else
        value = 0.0;
    end
elseif isfinite(E.time(index))
    value = abs(qesTime - E.time(index));
else
    value = 0.0;
end
end

function S = normalize_samplers(value)
if isstruct(value) && isfield(value,'qes') && numel(value) == 1
    value = value.qes;
end

if isstruct(value) && numel(value) == 1 && is_flat_struct(value)
    value = expand_flat_struct(value);
end

if istable(value)
    n = height(value);
    S.id = text_column(value,{'sampler_id','id'},n);
    S.type = text_column(value,{'sampler_type','type'},n);
    S.distance = numeric_column(value,{'distance_m','distance'},n);
    S.height = numeric_column(value,{'height_m','height','z'},n);
    S.azimuth = numeric_column(value,{'azimuth_deg','azimuth'},n);
    S.time = numeric_column(value,{'time_s','time'},n);
    S.windowStart = numeric_column(value, ...
        {'window_start','windowStart'},n);
    S.windowEnd = numeric_column(value,{'window_end','windowEnd'},n);
    S.sampleComplete = numeric_column(value,{'sample_complete'},n);
elseif isstruct(value)
    n = numel(value);
    S.id = text_field(value,{'sampler_id','id'},n);
    S.type = text_field(value,{'sampler_type','type'},n);
    S.distance = numeric_field(value,{'distance_m','distance'},n);
    S.height = numeric_field(value,{'height_m','height','z'},n);
    S.azimuth = numeric_field(value,{'azimuth_deg','azimuth'},n);
    S.time = numeric_field(value,{'time_s','time'},n);
    S.windowStart = numeric_field(value, ...
        {'window_start','windowStart'},n);
    S.windowEnd = numeric_field(value,{'window_end','windowEnd'},n);
    S.sampleComplete = numeric_field(value,{'sample_complete'},n);
else
    error('Sampler input must be a table or struct array.');
end

if isempty(S.id)
    S.id = cell(n,1);
    for k = 1:n
        S.id{k} = '';
    end
end
if isempty(S.type)
    S.type = cell(n,1);
    for k = 1:n
        S.type{k} = '';
    end
end
S.distance = fill_missing(S.distance,n);
S.height = fill_missing(S.height,n);
S.azimuth = fill_missing(S.azimuth,n);
S.time = fill_missing(S.time,n);
S.windowStart = fill_missing(S.windowStart,n);
S.windowEnd = fill_missing(S.windowEnd,n);
S.sampleComplete = fill_missing(S.sampleComplete,n);
end

function tf = is_flat_struct(value)
% A scalar struct whose sampler metadata are vector fields, as returned by
% read_qes_psd_sampler(), must be expanded before generic struct-array reads.
tf = flat_struct_count(value) > 1;
end

function expanded = expand_flat_struct(value)
fields = fieldnames(value);
n = flat_struct_count(value);

expanded = struct([]);
for i = 1:n
    entry = struct();
    for k = 1:numel(fields)
        fieldName = fields{k};
        fieldValue = value.(fieldName);
        if ischar(fieldValue)
            fieldEntry = fieldValue;
        elseif iscell(fieldValue)
            if numel(fieldValue) == n
                fieldEntry = fieldValue{i};
            elseif numel(fieldValue) == 1
                fieldEntry = fieldValue{1};
            else
                fieldEntry = fieldValue(min(i,numel(fieldValue)));
            end
        else
            if numel(fieldValue) == n
                fieldEntry = fieldValue(i);
            elseif numel(fieldValue) == 1
                fieldEntry = fieldValue(1);
            else
                fieldEntry = fieldValue(min(i,numel(fieldValue)));
            end
        end
        entry.(fieldName) = fieldEntry;
    end
    expanded(i,1) = entry; %#ok<AGROW>
end
end

function n = flat_struct_count(value)
names = {'distance_m','distance','height_m','height','z','time_s','time', ...
    'sampler_id','id','azimuth_deg','azimuth','sampling_time_s', ...
    'window_start','windowStart','window_end','windowEnd'};
n = 0;
for k = 1:numel(names)
    if isfield(value,names{k})
        fieldValue = value.(names{k});
        if ~ischar(fieldValue)
            n = max(n,numel(fieldValue));
        end
    end
end
if n == 0
    n = 1;
end
end

function values = numeric_column(T,names,n)
values = getcol(T,names,false);
if isempty(values)
    values = NaN(n,1);
else
    values = to_numeric(values(:));
end
end

function values = text_column(T,names,n)
values = getcol(T,names,false);
if isempty(values)
    values = cell(n,1);
    for k = 1:n
        values{k} = '';
    end
else
    values = normalize_text_cell(to_cellstr_compat(values));
    values = values(:);
end
end

function values = numeric_field(S,names,n)
values = [];
for k = 1:numel(names)
    if isfield(S,names{k})
        values = field_values(S,names{k},n);
        break;
    end
end
values = fill_missing(values,n);
end

function values = text_field(S,names,n)
values = [];
for k = 1:numel(names)
    if isfield(S,names{k})
        values = field_values(S,names{k},n);
        break;
    end
end
if isempty(values)
    values = cell(n,1);
    for k = 1:n
        values{k} = '';
    end
else
    values = normalize_text_cell(to_cellstr_compat(values));
    values = values(:);
end
end

function values = field_values(S,name,n)
values = cell(n,1);
for k = 1:n
    values{k} = S(k).(name);
end
end

function values = fill_missing(values,n)
if isempty(values)
    values = NaN(n,1);
    return;
end
values = to_numeric(values(:));
if numel(values) < n
    values(end+1:n,1) = NaN;
end
end

function values = to_numeric(values)
if iscell(values)
    converted = NaN(size(values));
    for k = 1:numel(values)
        item = values{k};
        if (isnumeric(item) || islogical(item)) && isscalar(item)
            converted(k) = double(item);
        elseif ischar(item) || isstring(item)
            converted(k) = str2double(item);
        end
    end
    values = converted;
end
values = double(values);
end

function values = normalize_text_cell(values)
for k = 1:numel(values)
    if ~ischar(values{k})
        values{k} = num2str(values{k});
    end
end
end
