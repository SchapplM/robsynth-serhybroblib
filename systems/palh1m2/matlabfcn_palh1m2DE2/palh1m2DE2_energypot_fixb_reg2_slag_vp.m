% Calculate inertial parameters regressor of potential energy for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh1m2DE2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energypot_fixb_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:06:22
% EndTime: 2020-05-02 21:06:22
% DurationCPUTime: 0.26s
% Computational Cost: add. (420->80), mult. (338->100), div. (0->0), fcn. (437->39), ass. (0->61)
t281 = sin(qJ(1));
t297 = g(2) * t281;
t287 = cos(qJ(1));
t298 = g(1) * t287;
t263 = t297 + t298;
t302 = g(3) * pkin(13);
t275 = pkin(18) - pkin(22);
t293 = -pkin(21) + t275;
t266 = -qJ(2) - qJ(3) + t293;
t274 = -qJ(2) + t275;
t253 = pkin(21) - atan2(cos(t274), -sin(t274));
t294 = -qJ(2) + t253;
t243 = -atan2(-sin(t266), cos(t266)) + t294;
t242 = -qJ(1) + t243;
t301 = sin(t242) / 0.2e1;
t241 = qJ(1) + t243;
t300 = -cos(t241) / 0.2e1;
t285 = cos(qJ(3));
t299 = pkin(5) * t285;
t296 = g(3) * t285;
t286 = cos(qJ(2));
t295 = t286 * t299 + pkin(15);
t276 = sin(pkin(19));
t277 = cos(pkin(19));
t279 = sin(qJ(3));
t256 = t285 * t276 + t279 * t277;
t257 = t279 * t276 - t285 * t277;
t246 = qJ(2) + atan2(t257, t256);
t280 = sin(qJ(2));
t247 = -g(3) * (t286 * pkin(1) + pkin(13)) + t263 * (t280 * pkin(1) - pkin(15));
t244 = sin(t246);
t245 = cos(t246);
t292 = -g(3) * t245 + t263 * t244;
t268 = -pkin(20) + t293;
t264 = sin(t268);
t265 = cos(t268);
t291 = g(3) * t264 + t263 * t265;
t290 = g(3) * t279 + t263 * t285;
t255 = -t263 * pkin(15) - t302;
t289 = cos(pkin(17));
t288 = cos(pkin(18));
t284 = cos(qJ(4));
t283 = sin(pkin(17));
t282 = sin(pkin(18));
t278 = sin(qJ(4));
t273 = cos(t275);
t272 = sin(t275);
t271 = pkin(5) * t279 + pkin(1);
t262 = -g(1) * t281 + g(2) * t287;
t259 = t283 * t282 + t289 * t288;
t258 = -t282 * t289 + t288 * t283;
t252 = -t279 * t263 + t296;
t249 = g(3) * t265 - t263 * t264;
t240 = cos(t242);
t237 = sin(t241);
t236 = g(3) * t258 + t263 * t259;
t235 = -g(3) * t259 + t263 * t258;
t234 = atan2(t257, -t256) + t246;
t233 = cos(t234);
t232 = sin(t234);
t1 = [0, 0, 0, 0, 0, 0, -t263, -t262, -g(3), -t302, 0, 0, 0, 0, 0, 0, -g(3) * t286 + t263 * t280, g(3) * t280 + t263 * t286, t262, t255, 0, 0, 0, 0, 0, 0, -t252 * t280 - t290 * t286, -t252 * t286 + t290 * t280, t262, t247, 0, 0, 0, 0, 0, 0, t291, t249, t262, (-pkin(5) * t296 + t263 * t271) * t280 - t295 * t298 + (-g(3) * t271 - t297 * t299) * t286 - pkin(15) * t297 - t302, 0, 0, 0, 0, 0, 0, t278 * t262 + t291 * t284, t284 * t262 - t291 * t278, -t249, (t263 * pkin(9) - g(3) * pkin(11)) * t265 + (g(3) * pkin(9) + t263 * pkin(11)) * t264 - (t286 * t271 + t280 * t299 + pkin(13)) * g(3) - t263 * (-t271 * t280 + t295), 0, 0, 0, 0, 0, 0, t235 * t286 + t280 * t236, -t280 * t235 + t236 * t286, t262, -(pkin(13) - pkin(16)) * g(3) + t263 * pkin(14), 0, 0, 0, 0, 0, 0, g(3) * t272 + t263 * t273, g(3) * t273 - t263 * t272, t262, t247, 0, 0, 0, 0, 0, 0, t292, g(3) * t244 + t263 * t245, t262, t255, 0, 0, 0, 0, 0, 0, g(3) * t233 - t263 * t232, -g(3) * t232 - t263 * t233, t262, t292 * pkin(2) + t255, 0, 0, 0, 0, 0, 0, cos(t243) * g(3) + (t240 / 0.2e1 + t300) * g(2) + (t301 + t237 / 0.2e1) * g(1), sin(t243) * g(3) + (t301 - t237 / 0.2e1) * g(2) + (-t240 / 0.2e1 + t300) * g(1), t262, (-sin(t294) * t263 + (-cos(t253) * t286 - sin(t253) * t280) * g(3)) * pkin(4) + t247;];
U_reg = t1;
