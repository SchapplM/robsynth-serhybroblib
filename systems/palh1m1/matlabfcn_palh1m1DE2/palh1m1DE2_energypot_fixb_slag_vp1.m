% Calculate potential energy for
% palh1m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1DE2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_energypot_fixb_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE2_energypot_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1DE2_energypot_fixb_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:02:14
% EndTime: 2020-04-14 20:02:20
% DurationCPUTime: 4.36s
% Computational Cost: add. (70637->195), mult. (105932->264), div. (5052->9), fcn. (67021->48), ass. (0->131)
t326 = pkin(12) + rSges(6,3);
t325 = -pkin(2) - pkin(13);
t324 = -pkin(2) + pkin(13);
t323 = -pkin(8) - pkin(3);
t322 = -pkin(8) + pkin(3);
t321 = -pkin(9) - pkin(11);
t320 = pkin(11) - pkin(9);
t274 = pkin(1) ^ 2;
t248 = sin(qJ(2));
t250 = sin(pkin(19));
t254 = cos(qJ(2));
t256 = cos(pkin(19));
t224 = t248 * t256 - t250 * t254;
t308 = pkin(7) * t224;
t294 = -0.2e1 * pkin(1) * t308 + t274;
t287 = pkin(7) ^ 2 + t294;
t291 = pkin(3) ^ 2 - pkin(8) ^ 2;
t212 = t287 + t291;
t217 = pkin(1) - t308;
t209 = sqrt(-((pkin(7) - t322) * (pkin(7) + t322) + t294) * ((pkin(7) - t323) * (pkin(7) + t323) + t294));
t225 = t248 * t250 + t254 * t256;
t303 = t209 * t225;
t205 = -pkin(7) * t303 + t212 * t217;
t319 = -t205 / 0.2e1;
t206 = pkin(7) * t212 * t225 + t209 * t217;
t318 = t206 / 0.2e1;
t317 = sin(pkin(23)) / 0.2e1;
t316 = sin(pkin(21)) / 0.2e1;
t253 = cos(qJ(3));
t315 = -t253 / 0.2e1;
t314 = cos(pkin(18)) / 0.2e1;
t313 = 0.1e1 / pkin(2) / 0.2e1;
t247 = sin(qJ(3));
t214 = 0.1e1 / t287;
t301 = t214 / pkin(3);
t202 = (t206 * t315 + t247 * t319) * t301;
t203 = (t205 * t315 + t247 * t318) * t301;
t238 = pkin(23) + pkin(22);
t231 = sin(t238);
t232 = cos(t238);
t184 = t202 * t232 + t203 * t231;
t311 = pkin(5) * t184;
t242 = sin(pkin(20));
t245 = cos(pkin(20));
t310 = pkin(6) * (-t242 * t253 - t245 * t247);
t309 = pkin(6) * (t242 * t247 - t245 * t253);
t296 = -0.2e1 * pkin(4) * t311 + pkin(5) ^ 2;
t288 = pkin(4) ^ 2 + t296;
t293 = pkin(9) ^ 2 - pkin(11) ^ 2;
t178 = t288 - t293;
t181 = -pkin(4) * t184 + pkin(5);
t176 = sqrt(-((pkin(4) - t320) * (pkin(4) + t320) + t296) * ((pkin(4) - t321) * (pkin(4) + t321) + t296));
t185 = -t202 * t231 + t203 * t232;
t305 = t176 * t185;
t174 = -pkin(4) * t305 + t178 * t181;
t175 = pkin(4) * t178 * t185 + t176 * t181;
t239 = qJ(2) + qJ(3);
t244 = cos(pkin(21));
t179 = 0.1e1 / t288;
t304 = t179 / pkin(11);
t170 = atan2((t174 * t316 + t175 * t244 / 0.2e1) * t304, (-t174 * t244 / 0.2e1 + t175 * t316) * t304) + t239;
t169 = cos(t170);
t307 = t169 * pkin(10);
t306 = t254 * pkin(1) + pkin(14);
t302 = t214 / pkin(8);
t246 = sin(qJ(4));
t249 = sin(qJ(1));
t300 = t246 * t249;
t255 = cos(qJ(1));
t299 = t246 * t255;
t252 = cos(qJ(4));
t298 = t249 * t252;
t297 = t252 * t255;
t243 = cos(pkin(23));
t189 = qJ(2) + atan2((t205 * t317 + t243 * t318) * t301, (t206 * t317 + t243 * t319) * t301);
t290 = pkin(1) * t310;
t216 = -0.2e1 * t290;
t267 = pkin(6) ^ 2;
t295 = t216 + t267;
t208 = sqrt(-((pkin(1) - t324) * (pkin(1) + t324) + t295) * ((pkin(1) - t325) * (pkin(1) + t325) + t295));
t272 = pkin(2) ^ 2;
t292 = t267 + t274;
t286 = -pkin(13) ^ 2 + t292;
t210 = t216 + t272 + t286;
t215 = -pkin(1) + t310;
t284 = 0.1e1 / (t216 + t292) * t313;
t200 = qJ(2) + atan2((-t208 * t215 + t210 * t309) * t284, (-t208 * t309 - t210 * t215) * t284);
t233 = sin(t239);
t289 = pkin(5) * t233 + t306;
t285 = t179 / pkin(9) / 0.2e1;
t283 = 0.1e1 / pkin(13) * t313;
t230 = -t248 * pkin(1) + pkin(16);
t186 = pkin(22) - t189;
t282 = -rSges(3,1) * t248 - rSges(3,2) * t254;
t234 = cos(t239);
t281 = rSges(4,1) * t234 - rSges(4,2) * t233;
t168 = sin(t170);
t280 = rSges(5,1) * t169 - rSges(5,2) * t168;
t187 = sin(t189);
t188 = cos(t189);
t279 = -rSges(8,1) * t187 - rSges(8,2) * t188;
t198 = sin(t200);
t199 = cos(t200);
t278 = -rSges(9,1) * t198 - rSges(9,2) * t199;
t211 = t287 - t291;
t218 = pkin(1) * t224 - pkin(7);
t204 = -pkin(1) * t303 - t211 * t218;
t207 = pkin(1) * t211 * t225 - t209 * t218;
t251 = sin(pkin(18));
t193 = atan2((t207 * t314 + t204 * t251 / 0.2e1) * t302, (t204 * t314 - t251 * t207 / 0.2e1) * t302);
t190 = sin(t193);
t191 = cos(t193);
t277 = rSges(7,1) * t191 - rSges(7,2) * t190 - pkin(15);
t196 = atan2(t208 * t283, (t272 - t286 + 0.2e1 * t290) * t283) + t200;
t194 = sin(t196);
t195 = cos(t196);
t276 = -pkin(2) * t198 + rSges(10,1) * t194 + rSges(10,2) * t195 + pkin(16);
t177 = t288 + t293;
t180 = -pkin(4) + t311;
t173 = -atan2((pkin(5) * t177 * t185 - t176 * t180) * t285, (-pkin(5) * t305 - t177 * t180) * t285) + t186;
t171 = sin(t173);
t172 = cos(t173);
t275 = -rSges(11,1) * t171 + rSges(11,2) * t172 + pkin(4) * sin(t186) + t230;
t236 = t255 * pkin(16);
t235 = t249 * pkin(16);
t228 = t255 * t230;
t227 = t249 * t230;
t226 = pkin(5) * t234 + t230;
t221 = t255 * t226;
t220 = t249 * t226;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t255 - rSges(2,2) * t249) + g(2) * (rSges(2,1) * t249 + rSges(2,2) * t255) + g(3) * (pkin(14) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t249 + t255 * t282 + t236) + g(2) * (-rSges(3,3) * t255 + t249 * t282 + t235) + g(3) * (rSges(3,1) * t254 - rSges(3,2) * t248 + pkin(14))) - m(4) * (g(1) * (rSges(4,3) * t249 + t255 * t281 + t228) + g(2) * (-rSges(4,3) * t255 + t249 * t281 + t227) + g(3) * (rSges(4,1) * t233 + rSges(4,2) * t234 + t306)) - m(5) * (g(1) * (rSges(5,3) * t249 + t255 * t280 + t221) + g(2) * (-rSges(5,3) * t255 + t249 * t280 + t220) + g(3) * (rSges(5,1) * t168 + rSges(5,2) * t169 + t289)) - m(6) * (g(1) * (t255 * t307 + t221 + (t169 * t297 + t300) * rSges(6,1) + (-t169 * t299 + t298) * rSges(6,2)) + g(2) * (t249 * t307 + t220 + (t169 * t298 - t299) * rSges(6,1) + (-t169 * t300 - t297) * rSges(6,2)) + g(3) * (-t326 * t169 + t289) + (g(3) * (rSges(6,1) * t252 - rSges(6,2) * t246 + pkin(10)) + (g(1) * t255 + g(2) * t249) * t326) * t168) - m(7) * (g(3) * (rSges(7,1) * t190 + rSges(7,2) * t191 + pkin(14) - pkin(17)) + (-g(2) * rSges(7,3) + g(1) * t277) * t255 + (g(1) * rSges(7,3) + g(2) * t277) * t249) - m(8) * (g(1) * (rSges(8,3) * t249 + t255 * t279 + t228) + g(2) * (-rSges(8,3) * t255 + t249 * t279 + t227) + g(3) * (rSges(8,1) * t188 - rSges(8,2) * t187 + t306)) - m(9) * (g(1) * (rSges(9,3) * t249 + t255 * t278 + t236) + g(2) * (-rSges(9,3) * t255 + t249 * t278 + t235) + g(3) * (rSges(9,1) * t199 - rSges(9,2) * t198 + pkin(14))) - m(10) * (g(3) * (pkin(2) * t199 - rSges(10,1) * t195 + rSges(10,2) * t194 + pkin(14)) + (-g(2) * rSges(10,3) + g(1) * t276) * t255 + (g(1) * rSges(10,3) + g(2) * t276) * t249) - m(11) * (g(3) * (pkin(4) * cos(t186) - t172 * rSges(11,1) - t171 * rSges(11,2) + t306) + (-g(2) * rSges(11,3) + g(1) * t275) * t255 + (g(1) * rSges(11,3) + g(2) * t275) * t249);
U = t1;
