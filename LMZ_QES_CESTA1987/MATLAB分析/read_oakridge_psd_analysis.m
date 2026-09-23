function PSD = read_oakridge_psd_analysis(fileName)
% Compatibility wrapper returning the PSD structure used by the
% aerosol-PSD evolution analysis.  The original experimental points are
% retained, and an additional 9-bin reference is produced by applying the
% same fixed QES binning rule used by AerosolPSDModel.

if nargin < 1
    PSD = read_oakridge_psd();
else
    PSD = read_oakridge_psd(fileName);
end

if ~isempty(PSD.exp.diameter_um)
    [PSD.exp.mass_fraction_9bin,PSD.exp.number_fraction_9bin] = ...
        map_to_qes_bins(PSD.exp.diameter_um,PSD.exp.mass_fraction, ...
        PSD.exp.number_fraction);
else
    PSD.exp.mass_fraction_9bin = zeros(1,9);
    PSD.exp.number_fraction_9bin = zeros(1,9);
end
end

function [mass9,number9] = map_to_qes_bins(diameter_um,mass,number)
% Fixed QES bin upper bounds in micrometres.  The final bin is open-ended.
upper_um = [0.4 0.7 1.1 2.1 3.3 4.7 5.8 9.0 Inf];
mass9 = zeros(1,9);
number9 = zeros(1,9);
m = mass(:);
n = number(:);
d = diameter_um(:);

for i = 1:numel(d)
    bin = find(d(i) < upper_um,1,'first');
    if isempty(bin)
        bin = 9;
    end
    if i <= numel(m)
        mass9(bin) = mass9(bin) + m(i);
    end
    if i <= numel(n)
        number9(bin) = number9(bin) + n(i);
    end
end

if sum(mass9) > 0
    mass9 = mass9 / sum(mass9);
end
if sum(number9) > 0
    number9 = number9 / sum(number9);
end
end
