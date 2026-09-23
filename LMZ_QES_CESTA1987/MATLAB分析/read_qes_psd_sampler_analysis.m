function PSD = read_qes_psd_sampler_analysis(fileName)
% Read aerosol_sampler.csv and expose one structure per QES row.

if nargin < 1
    flat = read_qes_psd_sampler();
else
    flat = read_qes_psd_sampler(fileName);
end

PSD = struct();
PSD.qes = struct('distance',{},'height',{},'azimuth',{},'time',{}, ...
    'sampler_id',{},'mass_fraction',{},'number_fraction',{}, ...
    'd50',{},'sigma_g',{},'windowStart',{},'windowEnd',{}, ...
    'sampling_time_s',{},'sample_complete',{},'dry_d_m',{},'wet_d_m',{}, ...
    'aerodynamic_d_m',{},'diameter_basis',{});

if ~isfield(flat,'qes') || isempty(flat.qes.time_s)
    return;
end

n = numel(flat.qes.time_s);
for i = 1:n
    PSD.qes(i).distance = flat.qes.distance_m(i);
    PSD.qes(i).height = flat.qes.height_m(i);
    if isfield(flat.qes,'azimuth_deg') && numel(flat.qes.azimuth_deg) >= i
        PSD.qes(i).azimuth = flat.qes.azimuth_deg(i);
    else
        PSD.qes(i).azimuth = NaN;
    end
    PSD.qes(i).time = flat.qes.time_s(i);
    PSD.qes(i).windowStart = flat.qes.window_start(i);
    PSD.qes(i).windowEnd = flat.qes.window_end(i);
    PSD.qes(i).sampling_time_s = flat.qes.sampling_time_s(i);
    PSD.qes(i).sample_complete = flat.qes.sample_complete(i);
    PSD.qes(i).dry_d_m = flat.qes.dry_d_m(i,:);
    PSD.qes(i).wet_d_m = flat.qes.wet_d_m(i,:);
    PSD.qes(i).aerodynamic_d_m = flat.qes.aerodynamic_d_m(i,:);
    PSD.qes(i).diameter_basis = flat.qes.diameter_basis;
    if isfield(flat.qes,'sampler_id') && numel(flat.qes.sampler_id) >= i
        PSD.qes(i).sampler_id = flat.qes.sampler_id{i};
    else
        PSD.qes(i).sampler_id = sprintf('qes_psd_%d',i);
    end
    PSD.qes(i).mass_fraction = flat.qes.mass_fraction(i,:);
    PSD.qes(i).number_fraction = flat.qes.number_fraction(i,:);
    PSD.qes(i).d50 = flat.qes.d50_um(i);
    PSD.qes(i).sigma_g = flat.qes.sigma_g(i);
end
end
