% Calculate inertial parameters regressor of potential energy for
% palh1m2DE1
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
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh1m2DE1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_energypot_fixb_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 21:02:18
% EndTime: 2020-05-01 21:02:18
% DurationCPUTime: 0.27s
% Computational Cost: add. (331->69), mult. (688->109), div. (0->0), fcn. (778->20), ass. (0->63)
t282 = sin(pkin(18));
t288 = cos(pkin(18));
t259 = -t282 * pkin(9) + t288 * pkin(11);
t260 = t288 * pkin(9) + t282 * pkin(11);
t270 = sin(pkin(22));
t274 = cos(pkin(22));
t247 = t259 * t274 + t270 * t260;
t271 = sin(pkin(21));
t275 = cos(pkin(21));
t290 = t270 * t259 - t260 * t274;
t304 = t247 * t275 - t290 * t271;
t303 = t247 * t271 + t290 * t275;
t302 = pkin(13) * g(3);
t285 = cos(qJ(3));
t301 = pkin(5) * t285;
t287 = cos(qJ(1));
t300 = g(1) * t287;
t281 = sin(qJ(1));
t299 = g(2) * t281;
t298 = t285 * g(3);
t286 = cos(qJ(2));
t297 = t286 * t301 + pkin(15);
t280 = sin(qJ(2));
t296 = pkin(1) * t280 - pkin(15);
t293 = g(3) * (t286 * pkin(1) + pkin(13));
t261 = t299 + t300;
t279 = sin(qJ(3));
t292 = t279 * g(3) + t261 * t285;
t273 = sin(pkin(19));
t277 = cos(pkin(19));
t249 = -t277 * g(3) + t261 * t273;
t250 = t273 * g(3) + t261 * t277;
t239 = t249 * t285 + t279 * t250;
t240 = t279 * t249 - t250 * t285;
t291 = t239 * t280 + t240 * t286;
t255 = t282 * t270 + t274 * t288;
t256 = t288 * t270 - t282 * t274;
t244 = t255 * t275 - t271 * t256;
t245 = t271 * t255 + t256 * t275;
t254 = -t261 * pkin(15) - t302;
t289 = cos(pkin(17));
t284 = cos(qJ(4));
t283 = sin(pkin(17));
t278 = sin(qJ(4));
t276 = cos(pkin(20));
t272 = sin(pkin(20));
t266 = pkin(5) * t279 + pkin(1);
t262 = -g(1) * t281 + g(2) * t287;
t258 = t283 * t282 + t289 * t288;
t257 = -t282 * t289 + t288 * t283;
t253 = -t279 * t261 + t298;
t252 = g(3) * t280 + t261 * t286;
t251 = -g(3) * t286 + t261 * t280;
t246 = t261 * t296 - t293;
t243 = -t280 * t253 - t292 * t286;
t242 = -t253 * t286 + t280 * t292;
t238 = g(3) * t257 + t261 * t258;
t237 = -g(3) * t258 + t261 * t257;
t236 = t244 * t272 + t245 * t276;
t235 = t244 * t276 - t272 * t245;
t233 = t235 * g(3) + t261 * t236;
t232 = -g(3) * t236 + t261 * t235;
t1 = [0, 0, 0, 0, 0, 0, -t261, -t262, -g(3), -t302, 0, 0, 0, 0, 0, 0, t251, t252, t262, t254, 0, 0, 0, 0, 0, 0, t243, t242, t262, t246, 0, 0, 0, 0, 0, 0, t232, t233, t262, (-pkin(5) * t298 + t261 * t266) * t280 - t297 * t300 + (-g(3) * t266 - t299 * t301) * t286 - pkin(15) * t299 - t302, 0, 0, 0, 0, 0, 0, t232 * t284 + t278 * t262, -t232 * t278 + t284 * t262, -t233, -(t286 * t266 - t303 * t272 + t304 * t276 + t280 * t301 + pkin(13)) * g(3) - t261 * (-t266 * t280 + t304 * t272 + t303 * t276 + t297), 0, 0, 0, 0, 0, 0, t237 * t286 + t238 * t280, -t237 * t280 + t238 * t286, t262, -(pkin(13) - pkin(16)) * g(3) + t261 * pkin(14), 0, 0, 0, 0, 0, 0, -g(3) * t256 + t261 * t255, g(3) * t255 + t261 * t256, t262, t246, 0, 0, 0, 0, 0, 0, t291, t239 * t286 - t280 * t240, t262, t254, 0, 0, 0, 0, 0, 0, t251, t252, t262, t291 * pkin(2) + t254, 0, 0, 0, 0, 0, 0, t243, t242, t262, -t293 + t261 * (t244 * pkin(4) + t296) - t245 * g(3) * pkin(4);];
U_reg = t1;
