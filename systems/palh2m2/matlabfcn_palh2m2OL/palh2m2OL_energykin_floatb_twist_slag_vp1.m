% Calculate kinetic energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
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
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m2OL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:18:28
% EndTime: 2020-05-03 01:18:31
% DurationCPUTime: 3.55s
% Computational Cost: add. (1908->352), mult. (1778->539), div. (0->0), fcn. (1572->16), ass. (0->192)
t234 = cos(qJ(3));
t209 = t234 * pkin(2) + pkin(4);
t230 = sin(qJ(2));
t235 = cos(qJ(2));
t229 = sin(qJ(3));
t287 = pkin(2) * t229;
t158 = t209 * t235 - t230 * t287 + pkin(1);
t312 = t158 * V_base(6);
t292 = t230 * t209;
t316 = t235 * t287 + t292;
t216 = V_base(6) + qJD(1);
t217 = V_base(4) * pkin(1);
t231 = sin(qJ(1));
t206 = t217 * t231;
t275 = t235 * V_base(4);
t271 = pkin(4) * t275;
t317 = t231 * t271 + t206;
t307 = t235 * pkin(4);
t210 = pkin(1) + t307;
t236 = cos(qJ(1));
t211 = V_base(5) * t236;
t315 = qJD(2) * t307 - t210 * t211 + V_base(3);
t205 = qJD(2) * t231 + V_base(4);
t181 = t235 * t229 + t230 * t234;
t283 = t316 * qJD(2);
t250 = -pkin(2) * t181 * qJD(3) - t283;
t314 = t250 * t236 + t316 * V_base(5) + V_base(1);
t268 = t235 * rSges(3,1) - t230 * rSges(3,2);
t251 = pkin(1) + t268;
t313 = t231 * rSges(3,3) + t251 * t236;
t310 = pkin(4) * t230;
t228 = sin(qJ(4));
t233 = cos(qJ(4));
t309 = pkin(5) * (t235 * (-t229 * t228 + t234 * t233) - t230 * (t234 * t228 + t229 * t233));
t308 = pkin(5) * ((t235 * t228 + t230 * t233) * t234 + t229 * (-t230 * t228 + t235 * t233));
t305 = t236 * rSges(3,3);
t304 = Icges(2,4) * t231;
t303 = Icges(3,4) * t230;
t302 = Icges(3,4) * t235;
t226 = qJ(2) + qJ(3);
t222 = sin(t226);
t301 = Icges(4,4) * t222;
t223 = cos(t226);
t300 = Icges(4,4) * t223;
t225 = qJ(4) + t226;
t212 = sin(t225);
t299 = Icges(5,4) * t212;
t213 = cos(t225);
t298 = Icges(5,4) * t213;
t214 = qJ(5) + t225;
t207 = sin(t214);
t297 = Icges(6,4) * t207;
t208 = cos(t214);
t296 = Icges(6,4) * t208;
t294 = t207 * t231;
t293 = t207 * t236;
t227 = sin(qJ(6));
t291 = t231 * t227;
t232 = cos(qJ(6));
t290 = t231 * t232;
t289 = t236 * t227;
t288 = t236 * t232;
t286 = pkin(3) * t208;
t220 = qJD(3) * t231;
t285 = qJD(6) * t207;
t284 = t158 * qJD(1);
t282 = qJD(2) + qJD(3);
t281 = t236 * t284 + V_base(2);
t279 = t216 * t309;
t273 = qJD(4) + t282;
t272 = t216 * t210;
t180 = t220 + t205;
t270 = t216 * t158;
t204 = -qJD(2) * t236 + V_base(5);
t267 = rSges(4,1) * t223 - rSges(4,2) * t222;
t266 = rSges(5,1) * t213 - rSges(5,2) * t212;
t265 = rSges(6,1) * t208 - rSges(6,2) * t207;
t201 = t232 * rSges(7,1) - rSges(7,2) * t227;
t264 = -rSges(7,3) * t207 + t201 * t208;
t165 = qJD(4) * t231 + t180;
t263 = Icges(3,1) * t235 - t303;
t262 = Icges(4,1) * t223 - t301;
t261 = Icges(5,1) * t213 - t299;
t260 = Icges(6,1) * t208 - t297;
t259 = -Icges(3,2) * t230 + t302;
t258 = -Icges(4,2) * t222 + t300;
t257 = -Icges(5,2) * t212 + t298;
t256 = -Icges(6,2) * t207 + t296;
t255 = Icges(3,5) * t235 - Icges(3,6) * t230;
t254 = Icges(4,5) * t223 - Icges(4,6) * t222;
t253 = Icges(5,5) * t213 - Icges(5,6) * t212;
t252 = Icges(6,5) * t208 - Icges(6,6) * t207;
t145 = qJD(5) * t231 + t165;
t173 = t231 * V_base(4) - t211 + t282;
t249 = pkin(2) * t173 * t223 + t315;
t144 = V_base(5) + (-qJD(5) - t273) * t236;
t248 = (-Icges(6,3) * t236 + t252 * t231) * t144 + (Icges(6,3) * t231 + t252 * t236) * t145 + (Icges(6,5) * t207 + Icges(6,6) * t208) * t216;
t164 = t273 * t236 - V_base(5);
t247 = (-Icges(5,3) * t236 + t253 * t231) * t164 - (Icges(5,3) * t231 + t253 * t236) * t165 - (Icges(5,5) * t212 + Icges(5,6) * t213) * t216;
t179 = -t282 * t236 + V_base(5);
t246 = (-Icges(4,3) * t236 + t254 * t231) * t179 + (Icges(4,3) * t231 + t254 * t236) * t180 + (Icges(4,5) * t222 + Icges(4,6) * t223) * t216;
t245 = (-Icges(3,3) * t236 + t255 * t231) * t204 + (Icges(3,3) * t231 + t255 * t236) * t205 + (Icges(3,5) * t230 + Icges(3,6) * t235) * t216;
t243 = -pkin(5) * (-qJD(4) - t173) * t213 + (t217 + t271) * t231 + t249;
t242 = -t164 * t308 + t314;
t119 = -Icges(6,6) * t236 + t256 * t231;
t120 = Icges(6,6) * t231 + t256 * t236;
t121 = -Icges(6,5) * t236 + t260 * t231;
t122 = Icges(6,5) * t231 + t260 * t236;
t160 = Icges(6,2) * t208 + t297;
t161 = Icges(6,1) * t207 + t296;
t241 = (-t120 * t207 + t122 * t208) * t145 + (-t119 * t207 + t121 * t208) * t144 + (-t160 * t207 + t161 * t208) * t216;
t128 = -Icges(5,6) * t236 + t257 * t231;
t129 = Icges(5,6) * t231 + t257 * t236;
t130 = -Icges(5,5) * t236 + t261 * t231;
t131 = Icges(5,5) * t231 + t261 * t236;
t168 = Icges(5,2) * t213 + t299;
t169 = Icges(5,1) * t212 + t298;
t240 = (-t129 * t212 + t131 * t213) * t165 - (-t128 * t212 + t130 * t213) * t164 + (-t168 * t212 + t169 * t213) * t216;
t138 = -Icges(4,6) * t236 + t258 * t231;
t139 = Icges(4,6) * t231 + t258 * t236;
t140 = -Icges(4,5) * t236 + t262 * t231;
t141 = Icges(4,5) * t231 + t262 * t236;
t176 = Icges(4,2) * t223 + t301;
t177 = Icges(4,1) * t222 + t300;
t239 = (-t139 * t222 + t141 * t223) * t180 + (-t138 * t222 + t140 * t223) * t179 + (-t176 * t222 + t177 * t223) * t216;
t150 = -Icges(3,6) * t236 + t259 * t231;
t151 = Icges(3,6) * t231 + t259 * t236;
t152 = -Icges(3,5) * t236 + t263 * t231;
t153 = Icges(3,5) * t231 + t263 * t236;
t192 = Icges(3,2) * t235 + t303;
t195 = Icges(3,1) * t230 + t302;
t238 = (-t151 * t230 + t153 * t235) * t205 + (-t150 * t230 + t152 * t235) * t204 + (-t192 * t230 + t195 * t235) * t216;
t237 = (-t181 * t220 - t229 * t275) * pkin(2) - t231 * t283 - t165 * t308 - V_base(4) * t292 + t281 + (t279 + t312) * t236;
t224 = Icges(2,4) * t236;
t202 = t236 * rSges(2,1) - t231 * rSges(2,2);
t200 = t227 * rSges(7,1) + t232 * rSges(7,2);
t199 = t231 * rSges(2,1) + t236 * rSges(2,2);
t198 = t230 * rSges(3,1) + t235 * rSges(3,2);
t197 = Icges(2,1) * t236 - t304;
t196 = Icges(2,1) * t231 + t224;
t194 = -Icges(2,2) * t231 + t224;
t193 = Icges(2,2) * t236 + t304;
t187 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t186 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t185 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t178 = t222 * rSges(4,1) + t223 * rSges(4,2);
t174 = qJD(6) * t208 + t216;
t170 = t212 * rSges(5,1) + t213 * rSges(5,2);
t162 = t207 * rSges(6,1) + t208 * rSges(6,2);
t157 = t208 * t288 - t291;
t156 = -t208 * t289 - t290;
t155 = t208 * t290 + t289;
t154 = -t208 * t291 + t288;
t147 = V_base(5) * rSges(2,3) - t216 * t199 + V_base(1);
t146 = -V_base(4) * rSges(2,3) + t216 * t202 + V_base(2);
t143 = t231 * rSges(4,3) + t267 * t236;
t142 = -t236 * rSges(4,3) + t267 * t231;
t133 = t231 * rSges(5,3) + t266 * t236;
t132 = -t236 * rSges(5,3) + t266 * t231;
t125 = t231 * rSges(6,3) + t265 * t236;
t124 = -t236 * rSges(6,3) + t265 * t231;
t123 = V_base(4) * t199 - V_base(5) * t202 + V_base(3);
t116 = t208 * rSges(7,3) + t201 * t207;
t113 = Icges(7,5) * t208 + (Icges(7,1) * t232 - Icges(7,4) * t227) * t207;
t112 = Icges(7,6) * t208 + (Icges(7,4) * t232 - Icges(7,2) * t227) * t207;
t111 = Icges(7,3) * t208 + (Icges(7,5) * t232 - Icges(7,6) * t227) * t207;
t110 = -t236 * t285 + t145;
t109 = -t231 * t285 + t144;
t107 = -t231 * t200 + t264 * t236;
t106 = t236 * t200 + t264 * t231;
t105 = t204 * t198 + V_base(1) + (-t251 * t231 + t305) * t216;
t104 = -t205 * t198 + t216 * t313 + V_base(2);
t103 = Icges(7,1) * t157 + Icges(7,4) * t156 - Icges(7,5) * t293;
t102 = Icges(7,1) * t155 + Icges(7,4) * t154 - Icges(7,5) * t294;
t101 = Icges(7,4) * t157 + Icges(7,2) * t156 - Icges(7,6) * t293;
t100 = Icges(7,4) * t155 + Icges(7,2) * t154 - Icges(7,6) * t294;
t99 = Icges(7,5) * t157 + Icges(7,6) * t156 - Icges(7,3) * t293;
t98 = Icges(7,5) * t155 + Icges(7,6) * t154 - Icges(7,3) * t294;
t97 = t268 * qJD(2) + V_base(3) + t206 + (t268 * t231 - t305) * V_base(4) - t313 * V_base(5);
t96 = t216 * t143 - t180 * t178 - t205 * t310 + t272 * t236 + V_base(2);
t95 = -t216 * t142 + t179 * t178 + t204 * t310 - t272 * t231 + V_base(1);
t94 = t180 * t142 - t179 * t143 + t315 + t317;
t93 = t165 * t132 + t164 * t133 + t249 + t317;
t92 = t145 * t124 - t144 * t125 + t243;
t91 = t216 * t133 - t165 * t170 + t250 * t231 + t236 * t312 - t316 * V_base(4) + t281;
t90 = -t216 * t132 - t164 * t170 + (-t284 - t312) * t231 + t314;
t89 = t216 * t125 - t145 * t162 + t237;
t88 = -t216 * t124 + t144 * t162 + (-t270 - t279) * t231 + t242;
t87 = t110 * t106 - t109 * t107 + (-t236 * t144 + t231 * t145) * t286 + t243;
t86 = t174 * t107 - t110 * t116 + (t208 * t236 * t216 - t207 * t145) * pkin(3) + t237;
t85 = t144 * t207 * pkin(3) - t174 * t106 + t109 * t116 + ((-t286 - t309) * t216 - t270) * t231 + t242;
t1 = m(1) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t92 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t91 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t110 * ((t156 * t101 + t157 * t103 - t99 * t293) * t110 + (t156 * t100 + t157 * t102 - t98 * t293) * t109 + (-t111 * t293 + t156 * t112 + t157 * t113) * t174) / 0.2e1 + t109 * ((t154 * t101 + t155 * t103 - t99 * t294) * t110 + (t154 * t100 + t155 * t102 - t98 * t294) * t109 + (-t111 * t294 + t154 * t112 + t155 * t113) * t174) / 0.2e1 + t205 * (t245 * t231 + t238 * t236) / 0.2e1 + t204 * (t238 * t231 - t245 * t236) / 0.2e1 + t180 * (t246 * t231 + t239 * t236) / 0.2e1 + t179 * (t239 * t231 - t246 * t236) / 0.2e1 - t164 * (t240 * t231 + t247 * t236) / 0.2e1 + t165 * (-t247 * t231 + t240 * t236) / 0.2e1 + t145 * (t248 * t231 + t241 * t236) / 0.2e1 + t144 * (t241 * t231 - t248 * t236) / 0.2e1 + m(2) * (t123 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(3) * (t104 ^ 2 + t105 ^ 2 + t97 ^ 2) / 0.2e1 + t174 * ((t98 * t109 + t99 * t110 + t111 * t174) * t208 + ((-t101 * t227 + t103 * t232) * t110 + (-t100 * t227 + t102 * t232) * t109 + (-t112 * t227 + t113 * t232) * t174) * t207) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + ((-t231 * t193 + t236 * t196 + Icges(1,4)) * V_base(5) + (-t231 * t194 + t236 * t197 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t236 * t193 + t231 * t196 + Icges(1,2)) * V_base(5) + (t236 * t194 + t231 * t197 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t235 * t151 + t230 * t153) * t205 + (t235 * t150 + t230 * t152) * t204 + (t223 * t139 + t222 * t141) * t180 + (t223 * t138 + t222 * t140) * t179 + (t213 * t129 + t212 * t131) * t165 - (t213 * t128 + t212 * t130) * t164 + (t208 * t120 + t207 * t122) * t145 + (t208 * t119 + t207 * t121) * t144 + (t208 * t160 + t207 * t161 + t213 * t168 + t212 * t169 + t223 * t176 + t222 * t177 + t235 * t192 + t230 * t195 + Icges(2,3)) * t216) * t216 / 0.2e1 + t216 * V_base(4) * (Icges(2,5) * t236 - Icges(2,6) * t231) + V_base(5) * t216 * (Icges(2,5) * t231 + Icges(2,6) * t236) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
