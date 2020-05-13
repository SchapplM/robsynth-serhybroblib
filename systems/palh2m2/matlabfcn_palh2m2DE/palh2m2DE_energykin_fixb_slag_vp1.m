% Calculate kinetic energy for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m2DE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2DE_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:24
% EndTime: 2020-05-03 01:06:25
% DurationCPUTime: 0.94s
% Computational Cost: add. (226->113), mult. (506->204), div. (0->0), fcn. (388->8), ass. (0->74)
t258 = sin(qJ(1));
t300 = t258 ^ 2;
t262 = cos(qJ(1));
t299 = t262 ^ 2;
t256 = sin(qJ(3));
t260 = cos(qJ(3));
t239 = Icges(5,5) * t256 + Icges(5,6) * t260;
t302 = t262 * t258;
t288 = t258 * qJD(1);
t301 = t262 * qJD(1);
t296 = pkin(4) * qJD(2);
t295 = pkin(5) * qJD(3);
t257 = sin(qJ(2));
t294 = Icges(3,4) * t257;
t261 = cos(qJ(2));
t293 = Icges(3,4) * t261;
t290 = qJD(2) * t258;
t289 = qJD(2) * t262;
t287 = t261 * pkin(4) + pkin(1);
t286 = t257 * t296;
t285 = t256 * t295;
t250 = pkin(2) + t287;
t284 = t262 * t286;
t248 = t261 * rSges(3,1) - t257 * rSges(3,2);
t283 = t260 * rSges(5,1) - t256 * rSges(5,2);
t282 = Icges(3,1) * t261 - t294;
t280 = -Icges(3,2) * t257 + t293;
t278 = Icges(3,5) * t261 - Icges(3,6) * t257;
t277 = Icges(5,5) * t260 - Icges(5,6) * t256;
t228 = -Icges(3,6) * t262 + t258 * t280;
t232 = -Icges(3,5) * t262 + t258 * t282;
t274 = t228 * t257 - t232 * t261;
t229 = Icges(3,6) * t258 + t262 * t280;
t233 = Icges(3,5) * t258 + t262 * t282;
t273 = -t229 * t257 + t233 * t261;
t242 = Icges(3,2) * t261 + t294;
t244 = Icges(3,1) * t257 + t293;
t271 = -t242 * t257 + t244 * t261;
t270 = t260 * pkin(5) + t250;
t269 = t250 + t283;
t268 = pkin(1) + t248;
t267 = -qJD(3) * (t256 * rSges(5,1) + t260 * rSges(5,2)) - t286;
t266 = -t285 - t286;
t265 = rSges(6,3) * qJD(1) + t266;
t263 = qJD(2) ^ 2;
t259 = cos(qJ(4));
t255 = sin(qJ(4));
t254 = qJD(1) + qJD(4);
t252 = t261 * t296;
t251 = rSges(4,1) + t287;
t249 = t262 * rSges(2,1) - t258 * rSges(2,2);
t247 = t258 * rSges(2,1) + t262 * rSges(2,2);
t246 = t257 * rSges(3,1) + t261 * rSges(3,2);
t240 = Icges(3,5) * t257 + Icges(3,6) * t261;
t238 = t260 * t295 + t252;
t237 = t238 ^ 2;
t236 = rSges(6,1) + t270;
t235 = pkin(3) + t270;
t234 = qJD(3) * t283 + t252;
t225 = Icges(3,3) * t258 + t262 * t278;
t224 = -Icges(3,3) * t262 + t258 * t278;
t223 = Icges(5,3) * t258 + t262 * t277;
t222 = -Icges(5,3) * t262 + t258 * t277;
t221 = (t262 * rSges(4,3) - t251 * t258) * qJD(1) - t284;
t220 = (t258 * rSges(4,3) + t251 * t262) * qJD(1) - t258 * t286;
t219 = -t236 * t288 + t262 * t265;
t218 = t236 * t301 + t258 * t265;
t217 = -t246 * t290 + (t258 * rSges(3,3) + t262 * t268) * qJD(1);
t216 = -t246 * t289 + (t262 * rSges(3,3) - t258 * t268) * qJD(1);
t215 = t267 * t258 + (t258 * rSges(5,3) + t262 * t269) * qJD(1);
t214 = t267 * t262 + (t262 * rSges(5,3) - t258 * t269) * qJD(1);
t213 = (t235 * qJD(1) + t254 * (rSges(7,1) * t259 - rSges(7,2) * t255)) * t262 + ((-t255 * rSges(7,1) - t259 * rSges(7,2)) * t254 + t266) * t258;
t212 = -t235 * t288 - t284 - t262 * t285 - t254 * ((rSges(7,1) * t258 + t262 * rSges(7,2)) * t259 + t255 * (t262 * rSges(7,1) - rSges(7,2) * t258));
t1 = m(6) * (t218 ^ 2 + t219 ^ 2 + t237) / 0.2e1 + m(5) * (t214 ^ 2 + t215 ^ 2 + t234 ^ 2) / 0.2e1 + m(7) * (t212 ^ 2 + t213 ^ 2 + t237) / 0.2e1 + m(4) * (pkin(4) ^ 2 * t261 ^ 2 * t263 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + m(3) * (t248 ^ 2 * t263 + t216 ^ 2 + t217 ^ 2) / 0.2e1 + ((t258 * t240 + t262 * t271) * qJD(1) + (t300 * t225 + (t274 * t262 + (-t224 + t273) * t258) * t262) * qJD(2)) * t290 / 0.2e1 - ((-t262 * t240 + t258 * t271) * qJD(1) + (t299 * t224 + (t273 * t258 + (-t225 + t274) * t262) * t258) * qJD(2)) * t289 / 0.2e1 + qJD(3) * t258 * (t239 * t288 + (-t222 * t302 + t300 * t223) * qJD(3)) / 0.2e1 - qJD(3) * t262 * (-t239 * t301 + (t299 * t222 - t223 * t302) * qJD(3)) / 0.2e1 + t254 ^ 2 * Icges(7,3) / 0.2e1 + (((t261 * t229 + t257 * t233) * t258 - (t261 * t228 + t257 * t232) * t262) * qJD(2) + (Icges(5,2) * t260 ^ 2 + t261 * t242 + t257 * t244 + (Icges(5,1) * t256 + 0.2e1 * Icges(5,4) * t260) * t256) * qJD(1) + (t299 + t300) * qJD(3) * t239) * qJD(1) / 0.2e1 + (m(2) * (t247 ^ 2 + t249 ^ 2) + Icges(6,2) + Icges(2,3) + Icges(4,2)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
