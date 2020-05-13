% Calculate inertial parameters regressor of potential energy for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh3m2DE2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energypot_fixb_reg2_slag_vp: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:22:44
% EndTime: 2020-05-07 04:22:44
% DurationCPUTime: 0.42s
% Computational Cost: add. (1068->133), mult. (304->131), div. (0->0), fcn. (392->68), ass. (0->79)
t264 = sin(qJ(1));
t266 = cos(qJ(1));
t231 = g(1) * t266 + g(2) * t264;
t294 = pkin(8) * g(2);
t293 = pkin(10) * g(2);
t292 = pkin(11) * g(3);
t276 = pkin(15) + qJ(2);
t253 = pkin(18) + t276;
t234 = pkin(17) + qJ(3) + t253;
t229 = atan2(sin(t253), -cos(t253));
t227 = pkin(17) - t229;
t275 = -qJ(2) + t227;
t194 = -atan2(-sin(t234), -cos(t234)) + t275;
t192 = qJ(1) + t194;
t291 = sin(t192) / 0.2e1;
t193 = -qJ(1) + t194;
t290 = cos(t193) / 0.2e1;
t232 = pkin(16) + t234;
t257 = qJ(2) + qJ(3);
t218 = atan2(-sin(t232), cos(t232)) + t257;
t216 = qJ(1) + t218;
t289 = sin(t216) / 0.2e1;
t217 = -qJ(1) + t218;
t205 = sin(t217);
t288 = t205 / 0.2e1;
t206 = cos(t216);
t287 = t206 / 0.2e1;
t286 = cos(t217) / 0.2e1;
t228 = qJ(2) + t229;
t226 = -qJ(1) + t228;
t285 = sin(t226) / 0.2e1;
t225 = qJ(1) + t228;
t284 = -cos(t225) / 0.2e1;
t254 = pkin(14) - t276;
t246 = qJ(1) + t254;
t283 = -sin(t246) / 0.2e1;
t247 = -qJ(1) + t254;
t282 = -cos(t247) / 0.2e1;
t256 = qJ(1) - t257;
t281 = sin(t256) / 0.2e1;
t255 = qJ(1) + t257;
t280 = cos(t255) / 0.2e1;
t279 = pkin(4) * sin(qJ(3));
t215 = -qJ(4) + t218;
t214 = qJ(4) + t218;
t208 = qJ(1) + t214;
t211 = -qJ(1) + t215;
t258 = qJ(1) + qJ(4);
t274 = -sin(t208) / 0.4e1 + sin(t211) / 0.4e1 - sin(t258) / 0.2e1;
t209 = -qJ(1) + t214;
t210 = qJ(1) + t215;
t259 = qJ(1) - qJ(4);
t273 = -sin(t209) / 0.4e1 + sin(t210) / 0.4e1 - sin(t259) / 0.2e1;
t272 = cos(t208) / 0.4e1 + cos(t211) / 0.4e1 + cos(t258) / 0.2e1;
t271 = cos(t209) / 0.4e1 + cos(t210) / 0.4e1 - cos(t259) / 0.2e1;
t270 = t231 * pkin(12);
t263 = sin(qJ(2));
t239 = pkin(1) * t263 + pkin(11);
t265 = cos(qJ(2));
t195 = -g(3) * t239 - t231 * (pkin(1) * t265 + pkin(12));
t241 = sin(t255);
t244 = cos(t256);
t269 = g(3) * sin(t257) + (t281 + t241 / 0.2e1) * g(2) + (t280 + t244 / 0.2e1) * g(1);
t268 = pkin(8) * g(1);
t267 = pkin(10) * g(1);
t261 = qJ(1) - qJ(2);
t260 = qJ(1) + qJ(2);
t245 = cos(qJ(3)) * pkin(4) - pkin(1);
t237 = cos(t246);
t236 = sin(t247);
t230 = -g(1) * t264 + g(2) * t266;
t224 = cos(t226);
t221 = sin(t225);
t213 = cos(t218);
t212 = sin(t218);
t190 = cos(t192);
t189 = sin(t193);
t187 = -g(3) * t213 + (t286 - t206 / 0.2e1) * g(2) + (t288 + t289) * g(1);
t1 = [0, 0, 0, 0, 0, 0, -t231, -t230, -g(3), -t292, 0, 0, 0, 0, 0, 0, -g(3) * t263 - t231 * t265, -g(3) * t265 + t231 * t263, t230, -t270 - t292, 0, 0, 0, 0, 0, 0, t269, g(3) * cos(t257) + (t280 - t244 / 0.2e1) * g(2) + (t281 - t241 / 0.2e1) * g(1), t230, t195, 0, 0, 0, 0, 0, 0, g(3) * t212 + (-t205 / 0.2e1 + t289) * g(2) + (t286 + t287) * g(1), -t187, t230, (g(3) * t279 + t231 * t245) * t265 + (t245 * t263 - pkin(11)) * g(3) - t231 * (t263 * t279 + pkin(12)), 0, 0, 0, 0, 0, 0, (sin(t214) / 0.2e1 + sin(t215) / 0.2e1) * g(3) + (t273 - t274) * g(2) + (t271 + t272) * g(1), (-cos(t215) / 0.2e1 + cos(t214) / 0.2e1) * g(3) + (-t271 + t272) * g(2) + (t273 + t274) * g(1), t187, (t268 + t293) * t286 + (t267 - t294) * t288 + (t268 - t293) * t287 + (t267 + t294) * t289 - t270 + (pkin(8) * t212 - pkin(10) * t213 - t239) * g(3) + ((-sin(t261) / 0.2e1 - sin(t260) / 0.2e1) * g(2) + (-cos(t261) / 0.2e1 - cos(t260) / 0.2e1) * g(1)) * pkin(1) + t269 * pkin(4), 0, 0, 0, 0, 0, 0, g(3) * sin(t254) + (t283 + t236 / 0.2e1) * g(2) + (-t237 / 0.2e1 + t282) * g(1), -g(3) * cos(t254) + (t237 / 0.2e1 + t282) * g(2) + (t283 - t236 / 0.2e1) * g(1), t230, -(pkin(11) + pkin(13)) * g(3) + t231 * pkin(6), 0, 0, 0, 0, 0, 0, -sin(t228) * g(3) + (-t221 / 0.2e1 + t285) * g(2) + (t284 - t224 / 0.2e1) * g(1), -cos(t228) * g(3) + (t284 + t224 / 0.2e1) * g(2) + (t221 / 0.2e1 + t285) * g(1), t230, t195, 0, 0, 0, 0, 0, 0, -sin(t194) * g(3) + (t291 - t189 / 0.2e1) * g(2) + (t190 / 0.2e1 + t290) * g(1), cos(t194) * g(3) + (-t190 / 0.2e1 + t290) * g(2) + (t291 + t189 / 0.2e1) * g(1), t230, (-cos(t275) * t231 + (-cos(t227) * t263 + sin(t227) * t265) * g(3)) * pkin(3) + t195;];
U_reg = t1;
