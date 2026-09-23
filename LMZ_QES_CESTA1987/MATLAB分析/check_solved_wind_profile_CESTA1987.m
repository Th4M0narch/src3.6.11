function out = check_solved_wind_profile_CESTA1987()
%CHECK_SOLVED_WIND_PROFILE_CESTA1987
% Verify the SOLVED QES-Winds field against the CESTA 1987 target profile
% at the release location.
%
% Expected location of this file:
%   I:\LMZ\QES\LMZ_QES_CESTA1987\MATLAB分析
%
% Reads:
%   ..\Output\lmz_windsOut.nc
%
% Checks:
%   1) nearest QES grid point to the source (550,100)
%   2) solved horizontal wind speed at z = 2,10,18,28 m
%   3) speed error relative to CESTA target values
%   4) solved flow-to azimuth and wind-from azimuth
%   5) near-ground profile at 0.25,0.5,1,2,3,10,18,28 m
%   6) horizontal-domain statistics at the four CESTA heights
%
% Direction convention used here:
%   x = east, y = north
%   flow-to azimuth = atan2(uEast,vNorth), clockwise from north
%   wind-from azimuth = flow-to + 180 deg
%
% This function does not modify any QES file.

cfg = CESTA1987_config();
out = struct();

ncFile = cfg.windsVizNc;
if ~exist(ncFile,'file')
    error('Wind NetCDF not found: %s',ncFile);
end

fprintf('\n------------------------------------------------------------\n');
fprintf('Solved QES-Winds profile check at CESTA source\n');
fprintf('NetCDF : %s\n',ncFile);
fprintf('Target source XY : %.3f, %.3f m\n',cfg.sourceXYZ(1),cfg.sourceXYZ(2));
fprintf('------------------------------------------------------------\n');

% -------- Read coordinates and fields --------
x = double(ncread(ncFile,'x')); x=x(:);
y = double(ncread(ncFile,'y')); y=y(:);
z = double(ncread(ncFile,'z')); z=z(:);

u = double(ncread(ncFile,'u'));
v = double(ncread(ncFile,'v'));

% Drop singleton time dimension if present.
u = squeeze(u);
v = squeeze(v);

if ndims(u) ~= 3 || ndims(v) ~= 3
    error('Expected u/v to reduce to 3-D arrays after squeeze. Got u=%s v=%s',...
        mat2str(size(u)),mat2str(size(v)));
end

% NetCDF inventory supplied for this case is [x y z time].
if size(u,1) ~= numel(x) || size(u,2) ~= numel(y) || size(u,3) ~= numel(z)
    error(['u dimensions do not match x/y/z as [x y z]. ' ...
           'u=%s, nx=%d ny=%d nz=%d'],...
           mat2str(size(u)),numel(x),numel(y),numel(z));
end

[~,ix] = min(abs(x-cfg.sourceXYZ(1)));
[~,iy] = min(abs(y-cfg.sourceXYZ(2)));

fprintf('Nearest grid XY   : x=%.6f m (ix=%d), y=%.6f m (iy=%d)\n',...
    x(ix),ix,y(iy),iy);
fprintf('XY offsets        : dx=%+.6f m, dy=%+.6f m\n',...
    x(ix)-cfg.sourceXYZ(1),y(iy)-cfg.sourceXYZ(2));

% Extract full vertical profile at the source cell.
us = squeeze(u(ix,iy,:));
vs = squeeze(v(ix,iy,:));
speed = hypot(us,vs);

flowTo = mod(atan2d(us,vs),360);
windFrom = mod(flowTo+180,360);

% -------- Four CESTA target heights --------
zt = cfg.windZ_m(:);
Ut = cfg.windSpeed_m_s(:);
nT = numel(zt);

zGrid = nan(nT,1);
Uq = nan(nT,1);
uQ = nan(nT,1);
vQ = nan(nT,1);
flowQ = nan(nT,1);
fromQ = nan(nT,1);
speedErr = nan(nT,1);
speedErrPct = nan(nT,1);
flowErrDeg = nan(nT,1);
fromErrDeg = nan(nT,1);
domainMean = nan(nT,1);
domainStd = nan(nT,1);
domainMin = nan(nT,1);
domainMax = nan(nT,1);

for k=1:nT
    [~,iz] = min(abs(z-zt(k)));
    zGrid(k)=z(iz);
    Uq(k)=speed(iz);
    uQ(k)=us(iz);
    vQ(k)=vs(iz);
    flowQ(k)=flowTo(iz);
    fromQ(k)=windFrom(iz);
    speedErr(k)=Uq(k)-Ut(k);
    speedErrPct(k)=100*speedErr(k)/Ut(k);
    flowErrDeg(k)=angleDiff(flowQ(k),cfg.plumeAzimuth_deg);
    fromErrDeg(k)=angleDiff(fromQ(k),cfg.windFrom_deg);

    sl = hypot(u(:,:,iz),v(:,:,iz));
    ss = sl(isfinite(sl));
    domainMean(k)=mean(ss);
    domainStd(k)=std(ss);
    domainMin(k)=min(ss);
    domainMax(k)=max(ss);
end

T = table(zt,zGrid,Ut,Uq,speedErr,speedErrPct,uQ,vQ,...
          flowQ,flowErrDeg,fromQ,fromErrDeg,...
          domainMean,domainStd,domainMin,domainMax,...
    'VariableNames',{'z_target_m','z_grid_m','U_target_m_s','U_QES_m_s',...
    'speed_error_m_s','speed_error_pct','u_m_s','v_m_s',...
    'flow_to_deg','flow_to_error_deg','wind_from_deg','wind_from_error_deg',...
    'slice_mean_U_m_s','slice_std_U_m_s','slice_min_U_m_s','slice_max_U_m_s'});

disp(T);

outFile = fullfile(cfg.resultsDir,'solved_wind_profile_check.csv');
writetable(T,outFile);
fprintf('Saved four-height check: %s\n',outFile);

% Summary indicators.
out.maxAbsSpeedErrorPct = max(abs(speedErrPct));
out.meanAbsSpeedErrorPct = mean(abs(speedErrPct));
out.maxAbsFlowToErrorDeg = max(abs(flowErrDeg));
out.maxAbsWindFromErrorDeg = max(abs(fromErrDeg));

fprintf('\nFour-height summary:\n');
fprintf('  max |speed error|     = %.4f %%\n',out.maxAbsSpeedErrorPct);
fprintf('  mean |speed error|    = %.4f %%\n',out.meanAbsSpeedErrorPct);
fprintf('  max |flow-to error|   = %.4f deg\n',out.maxAbsFlowToErrorDeg);
fprintf('  max |wind-from error| = %.4f deg\n',out.maxAbsWindFromErrorDeg);

% Informational criteria: not a model "pass/fail" standard from CESTA.
% These are convenient engineering flags only.
speedFlag = out.maxAbsSpeedErrorPct <= 10;
dirFlag   = out.maxAbsFlowToErrorDeg <= 5;
fprintf('  engineering speed flag (<=10%%): %s\n',tf(speedFlag));
fprintf('  engineering dir flag   (<=5deg): %s\n',tf(dirFlag));
fprintf('  NOTE: these two thresholds are analysis flags, not CESTA acceptance criteria.\n');

% -------- Near-ground and extended profile --------
zCheck = [0.25 0.5 1 2 3 10 18 28]';
nC=numel(zCheck);
zActual=nan(nC,1); UCheck=nan(nC,1); uCheck=nan(nC,1); vCheck=nan(nC,1);
flowCheck=nan(nC,1); fromCheck=nan(nC,1);
for k=1:nC
    [~,iz]=min(abs(z-zCheck(k)));
    zActual(k)=z(iz);
    UCheck(k)=speed(iz);
    uCheck(k)=us(iz);
    vCheck(k)=vs(iz);
    flowCheck(k)=flowTo(iz);
    fromCheck(k)=windFrom(iz);
end
Tnear=table(zCheck,zActual,UCheck,uCheck,vCheck,flowCheck,fromCheck,...
    'VariableNames',{'z_requested_m','z_grid_m','U_QES_m_s','u_m_s','v_m_s',...
                     'flow_to_deg','wind_from_deg'});
writetable(Tnear,fullfile(cfg.resultsDir,'solved_wind_profile_near_ground.csv'));

fprintf('\nNear-ground/source profile:\n');
disp(Tnear);

% -------- Plot 1: solved profile vs target --------
fig=figure('Name','Solved CESTA wind profile check','Color','w');
plot(speed,z,'-','LineWidth',1.2); hold on;
plot(Ut,zt,'ko','MarkerSize',7,'LineWidth',1.2);
plot(Uq,zGrid,'x','MarkerSize',8,'LineWidth',1.2);
grid on;
xlabel('Horizontal wind speed (m/s)');
ylabel('Height z (m)');
title('CESTA 1987: solved QES-Winds profile at release location');
legend('Solved QES profile','CESTA target','Nearest-grid comparison','Location','best');
ylim([0 max(30,max(z))]);
savefig_png(fig,fullfile(cfg.resultsDir,'solved_wind_profile_source.png'));

% -------- Plot 2: percent speed error --------
fig=figure('Name','CESTA wind speed error','Color','w');
plot(speedErrPct,zt,'o-','LineWidth',1.2);
xline_compat(0);
grid on;
xlabel('Speed error (%)');
ylabel('Height z (m)');
title('QES-Winds speed error relative to CESTA target');
savefig_png(fig,fullfile(cfg.resultsDir,'solved_wind_profile_speed_error.png'));

% -------- Plot 3: directions --------
fig=figure('Name','CESTA wind direction check','Color','w');
plot(flowQ,zt,'o-','LineWidth',1.2); hold on;
plot(repmat(cfg.plumeAzimuth_deg,size(zt)),zt,'k--','LineWidth',1.0);
plot(fromQ,zt,'s-','LineWidth',1.2);
plot(repmat(cfg.windFrom_deg,size(zt)),zt,'k:','LineWidth',1.0);
grid on;
xlabel('Azimuth (deg clockwise from north)');
ylabel('Height z (m)');
title('Solved wind direction at CESTA source');
legend('QES flow-to','Target flow-to','QES wind-from','Target wind-from','Location','best');
savefig_png(fig,fullfile(cfg.resultsDir,'solved_wind_direction_source.png'));

% Save whole source profile for later D3 documentation.
Tall=table(z,speed,us,vs,flowTo,windFrom,...
    'VariableNames',{'z_m','U_horizontal_m_s','u_m_s','v_m_s','flow_to_deg','wind_from_deg'});
writetable(Tall,fullfile(cfg.resultsDir,'solved_wind_profile_full_source.csv'));

out.sourceGridXY=[x(ix) y(iy)];
out.targetTable=T;
out.nearGroundTable=Tnear;
out.fullProfile=Tall;
out.speedEngineeringFlag=speedFlag;
out.directionEngineeringFlag=dirFlag;

fprintf('Solved wind-profile check completed.\n');
fprintf('------------------------------------------------------------\n');
end

function d = angleDiff(a,b)
% Signed shortest angular difference a-b in degrees, in [-180,180).
d = mod(a-b+180,360)-180;
end

function s=tf(x)
if x, s='PASS'; else, s='CHECK'; end
end

function xline_compat(x)
% Older-MATLAB-compatible vertical reference line.
yl=ylim;
plot([x x],yl,'k--','LineWidth',0.8);
ylim(yl);
end
