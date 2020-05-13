% Calculate kinetic energy for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m1DE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1DE_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:00
% EndTime: 2020-05-02 23:52:00
% DurationCPUTime: 0.72s
% Computational Cost: add. (322->126), mult. (531->223), div. (0->0), fcn. (431->10), ass. (0->82)
t261 = cos(qJ(3));
t250 = t261 * pkin(3) + pkin(2);
t257 = sin(qJ(3));
t258 = sin(qJ(2));
t262 = cos(qJ(2));
t279 = -t258 * t257 * pkin(3) + t250 * t262 + pkin(1);
t298 = rSges(5,1) + t279;
t295 = pkin(2) * qJD(2);
t294 = Icges(3,4) * t258;
t293 = Icges(3,4) * t262;
t255 = qJ(2) + qJ(3);
t252 = sin(t255);
t292 = Icges(4,4) * t252;
t253 = cos(t255);
t291 = Icges(4,4) * t253;
t290 = t262 * t257;
t259 = sin(qJ(1));
t289 = qJD(2) * t259;
t263 = cos(qJ(1));
t288 = qJD(2) * t263;
t287 = (pkin(3) * t290 + t258 * t250) * qJD(2);
t286 = t259 * qJD(1);
t285 = qJD(2) + qJD(3);
t284 = t258 * t295;
t283 = pkin(3) * (t258 * t261 + t290) * qJD(3);
t282 = t262 * t295;
t247 = -t262 * rSges(3,1) + t258 * rSges(3,2);
t280 = rSges(4,1) * t253 - rSges(4,2) * t252;
t278 = Icges(3,1) * t262 - t294;
t277 = Icges(4,1) * t253 - t292;
t276 = -Icges(3,2) * t258 + t293;
t275 = -Icges(4,2) * t252 + t291;
t274 = Icges(3,5) * t262 - Icges(3,6) * t258;
t273 = Icges(4,5) * t253 - Icges(4,6) * t252;
t227 = Icges(3,6) * t263 + t276 * t259;
t229 = Icges(3,5) * t263 + t278 * t259;
t272 = -t227 * t258 + t229 * t262;
t228 = -Icges(3,6) * t259 + t276 * t263;
t230 = -Icges(3,5) * t259 + t278 * t263;
t271 = t228 * t258 - t230 * t262;
t243 = -Icges(3,2) * t262 - t294;
t244 = -Icges(3,1) * t258 - t293;
t270 = -t243 * t258 + t244 * t262;
t269 = pkin(1) - t247;
t268 = -t283 - t287;
t240 = t285 * t259;
t241 = t285 * t263;
t267 = (-Icges(4,5) * t252 - Icges(4,6) * t253) * qJD(1) + (Icges(4,3) * t263 + t273 * t259) * t241 - (-Icges(4,3) * t259 + t273 * t263) * t240;
t266 = -rSges(5,3) * qJD(1) + t268;
t219 = Icges(4,6) * t263 + t275 * t259;
t220 = -Icges(4,6) * t259 + t275 * t263;
t221 = Icges(4,5) * t263 + t277 * t259;
t222 = -Icges(4,5) * t259 + t277 * t263;
t236 = -Icges(4,2) * t253 - t292;
t237 = -Icges(4,1) * t252 - t291;
t265 = -(-t220 * t252 + t222 * t253) * t240 + (-t219 * t252 + t221 * t253) * t241 + (-t236 * t252 + t237 * t253) * qJD(1);
t260 = cos(qJ(4));
t256 = sin(qJ(4));
t254 = qJD(1) + qJD(4);
t251 = t262 * pkin(2) + pkin(1);
t248 = t263 * rSges(2,1) - t259 * rSges(2,2);
t246 = t259 * rSges(2,1) + t263 * rSges(2,2);
t245 = -t258 * rSges(3,1) - t262 * rSges(3,2);
t242 = -Icges(3,5) * t258 - Icges(3,6) * t262;
t238 = -t252 * rSges(4,1) - t253 * rSges(4,2);
t234 = -pkin(3) * t285 * t253 - t282;
t233 = t234 ^ 2;
t231 = pkin(4) + t279;
t226 = -Icges(3,3) * t259 + t274 * t263;
t225 = Icges(3,3) * t263 + t274 * t259;
t224 = -t259 * rSges(4,3) + t280 * t263;
t223 = t263 * rSges(4,3) + t280 * t259;
t216 = t245 * t289 + (-t259 * rSges(3,3) + t269 * t263) * qJD(1);
t215 = t245 * t288 + (-t263 * rSges(3,3) - t269 * t259) * qJD(1);
t214 = -t259 * t284 + t240 * t238 + (t251 * t263 + t224) * qJD(1);
t213 = -t263 * t284 + t241 * t238 + (-t251 * t259 - t223) * qJD(1);
t212 = t266 * t263 - t298 * t286;
t211 = t298 * t263 * qJD(1) + t266 * t259;
t210 = -t240 * t223 - t241 * t224 - t282;
t209 = (t231 * qJD(1) + t254 * (rSges(6,1) * t260 - rSges(6,2) * t256)) * t263 + (t254 * (-t256 * rSges(6,1) - t260 * rSges(6,2)) + t268) * t259;
t208 = -t231 * t286 - t263 * t287 - t263 * t283 - t254 * ((rSges(6,1) * t259 + t263 * rSges(6,2)) * t260 + t256 * (t263 * rSges(6,1) - t259 * rSges(6,2)));
t1 = m(3) * (t247 ^ 2 * qJD(2) ^ 2 + t215 ^ 2 + t216 ^ 2) / 0.2e1 - ((-t259 * t242 + t270 * t263) * qJD(1) + (t259 ^ 2 * t226 + (t272 * t263 + (-t225 + t271) * t259) * t263) * qJD(2)) * t289 / 0.2e1 + ((t263 * t242 + t270 * t259) * qJD(1) + (t263 ^ 2 * t225 + (t271 * t259 + (-t226 + t272) * t263) * t259) * qJD(2)) * t288 / 0.2e1 + m(4) * (t210 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 - t240 * (-t267 * t259 + t265 * t263) / 0.2e1 + t241 * (t265 * t259 + t267 * t263) / 0.2e1 + m(5) * (t211 ^ 2 + t212 ^ 2 + t233) / 0.2e1 + m(6) * (t208 ^ 2 + t209 ^ 2 + t233) / 0.2e1 + t254 ^ 2 * Icges(6,3) / 0.2e1 + ((-(-t262 * t228 - t258 * t230) * t259 + (-t262 * t227 - t258 * t229) * t263) * qJD(2) - (-t253 * t220 - t252 * t222) * t240 + (-t253 * t219 - t252 * t221) * t241 + (-t253 * t236 - t252 * t237 - t262 * t243 - t258 * t244) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t246 ^ 2 + t248 ^ 2) + Icges(2,3) + Icges(5,2)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
