% Jacobian time derivative of explicit kinematic constraints of
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% WD [16x4]
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)
% [ParkChoPlo1999] Park, FC and Choi, Jihyeon and Ploen, SR: Symbolic formulation of closed chain dynamics in independent coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function WD = palh1m2DE2_kinconstr_expl_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_kinconstr_expl_jacobianD_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_kinconstr_expl_jacobianD_mdh_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_kinconstr_expl_jacobianD_mdh_sym_varpar: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:25
% EndTime: 2020-05-02 20:57:27
% DurationCPUTime: 1.31s
% Computational Cost: add. (4388->78), mult. (8552->167), div. (280->24), fcn. (10944->16), ass. (0->107)
t254 = sin(qJ(2));
t257 = cos(qJ(3));
t258 = cos(qJ(2));
t306 = sin(qJ(3));
t236 = -t254 * t306 + t257 * t258;
t239 = t257 * t254 + t258 * t306;
t255 = sin(pkin(18));
t259 = cos(pkin(18));
t194 = t236 * t259 + t255 * t239;
t198 = t255 * t236 - t239 * t259;
t251 = sin(pkin(20));
t252 = cos(pkin(20));
t159 = t194 * t252 + t251 * t198;
t249 = pkin(22) + pkin(21);
t247 = sin(t249);
t248 = cos(t249);
t309 = t251 * t194 - t198 * t252;
t134 = t159 * t248 - t247 * t309;
t132 = 0.1e1 / t134 ^ 2;
t342 = t159 * t247 + t248 * t309;
t352 = t132 * t342 ^ 2;
t310 = t194 * t248 + t198 * t247;
t151 = 0.1e1 / t310 ^ 2;
t157 = -t194 * t247 + t198 * t248;
t348 = t151 * t157 ^ 2;
t131 = 0.1e1 / t134;
t312 = qJD(2) + qJD(3);
t203 = t312 * t236;
t204 = t312 * t239;
t172 = -t203 * t259 - t255 * t204;
t173 = t255 * t203 - t204 * t259;
t146 = t251 * t172 + t173 * t252;
t311 = -t172 * t252 + t251 * t173;
t343 = t146 * t248 - t247 * t311;
t351 = t343 * t131;
t350 = t351 * t352;
t347 = t342 * (t146 * t247 + t248 * t311);
t142 = t172 * t247 + t173 * t248;
t150 = 0.1e1 / t310;
t332 = t142 * t150;
t346 = t332 * t348;
t344 = t157 * (t172 * t248 - t173 * t247);
t256 = sin(pkin(17));
t260 = cos(pkin(17));
t238 = t255 * t260 - t259 * t256;
t241 = t256 * t255 + t260 * t259;
t314 = t238 * t258 - t241 * t254;
t192 = 0.1e1 / t314 ^ 2;
t341 = t192 * t314;
t200 = t254 * t238 + t241 * t258;
t339 = t192 * t200 ^ 2;
t250 = sin(pkin(22));
t305 = cos(pkin(22));
t231 = t255 * t250 + t305 * t259;
t232 = t259 * t250 - t255 * t305;
t186 = t254 * t231 + t232 * t258;
t182 = 0.1e1 / t186 ^ 2;
t340 = t182 * t186;
t315 = t231 * t258 - t232 * t254;
t328 = t182 * t315 ^ 2;
t191 = 0.1e1 / t314;
t338 = t191 * t339;
t253 = cos(pkin(19));
t304 = sin(pkin(19));
t227 = t306 * t253 + t257 * t304;
t308 = 0.1e1 / t227 ^ 2;
t337 = t227 * t308;
t336 = t200 * t341;
t229 = t257 * t253 - t306 * t304;
t226 = t229 ^ 2;
t288 = t226 * t308;
t317 = 0.1e1 + t288;
t333 = 0.1e1 / t317 ^ 2;
t205 = 0.1e1 / t317;
t329 = qJD(2) * t336;
t181 = 0.1e1 / t186;
t327 = t181 * t328;
t326 = t315 * t340;
t323 = qJD(3) * t337;
t223 = 0.1e1 / t227;
t316 = qJD(2) * t326;
t280 = t229 * t323;
t237 = -t254 * t259 + t258 * t255;
t233 = 0.1e1 / t237 ^ 2;
t313 = qJD(2) * t233;
t307 = t223 * t308;
t216 = t229 * qJD(3);
t292 = t216 * t307;
t300 = (-t226 * t292 - t280) * t333;
t240 = t254 * t255 + t258 * t259;
t235 = t240 ^ 2;
t213 = t235 * t233 + 0.1e1;
t290 = t240 / t237 * t313;
t291 = t237 * t313;
t299 = (t235 * t290 + t240 * t291) / t213 ^ 2;
t289 = t307 * t226;
t211 = 0.1e1 / t213;
t180 = t200 * qJD(2);
t176 = t315 * qJD(2);
t174 = 0.1e1 + t339;
t166 = 0.1e1 + t328;
t140 = 0.1e1 + t348;
t125 = 0.1e1 + t352;
t122 = 0.2e1 * (t223 * t227 + t288) * t333 * (-t216 * t289 - t280) + (0.2e1 * t280 + (-t223 + 0.2e1 * t289 + t337) * t216) * t205;
t119 = 0.2e1 * (-t150 * t310 - t348) / t140 ^ 2 * (t151 * t344 - t346) + (t332 - 0.2e1 * t346 + (-t310 * t142 + 0.2e1 * t344) * t151) / t140;
t118 = 0.2e1 * (t131 * t134 + t352) / t125 ^ 2 * (t132 * t347 - t350) + (-t351 + 0.2e1 * t350 + (t134 * t343 - 0.2e1 * t347) * t132) / t125;
t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, t118, t118, 0; 0, 0, 0, 0; 0, 0.2e1 * (-t191 * t314 - t339) / t174 ^ 2 * (t180 * t338 + t329) + (t329 - (-0.2e1 * t338 - t341) * t180 + (-t200 * t191 + t336) * qJD(2)) / t174, 0, 0; 0, 0.2e1 * (t181 * t186 + t328) / t166 ^ 2 * (-t176 * t327 - t316) + (t316 + (0.2e1 * t327 + t340) * t176 + (-t181 * t315 + t326) * qJD(2)) / t166, 0, 0; 0, 0, -0.2e1 * t300 + 0.2e1 * (-t205 * t323 + (-t205 * t292 - t300 * t308) * t229) * t229, 0; 0, 0, t122, 0; 0, t119, t119, 0; 0, -0.2e1 * t299 + 0.2e1 * (t211 * t291 + (t211 * t290 - t233 * t299) * t240) * t240, 0, 0; 0, 0, t122, 0; 0, t119, t119, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
WD = t1;
