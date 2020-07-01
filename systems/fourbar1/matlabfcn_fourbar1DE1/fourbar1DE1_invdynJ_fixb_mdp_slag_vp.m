% Calculate vector of inverse dynamics joint torques for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1DE1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:30
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1DE1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_invdynJ_fixb_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_invdynJ_fixb_mdp_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbar1DE1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_invdynJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbar1DE1_invdynJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:29:57
% EndTime: 2020-06-26 17:30:05
% DurationCPUTime: 1.70s
% Computational Cost: add. (1858->128), mult. (3019->229), div. (76->10), fcn. (584->4), ass. (0->106)
t213 = pkin(2) ^ 2;
t216 = pkin(1) ^ 2;
t201 = cos(qJ(1));
t277 = pkin(2) * t201;
t257 = pkin(1) * t277;
t262 = -0.2e1 * t257 + t216;
t182 = t213 + t262;
t180 = 0.1e1 / t182 ^ 2;
t272 = t180 * t213;
t207 = pkin(4) ^ 2;
t209 = pkin(3) ^ 2;
t189 = t207 + t209;
t202 = pkin(3) + pkin(4);
t282 = pkin(4) - pkin(3);
t269 = (t282 ^ 2 * t202 ^ 2);
t285 = 0.2e1 * t182 * t189 - (2 * t269);
t206 = qJD(1) ^ 2;
t200 = sin(qJ(1));
t280 = pkin(1) * t200;
t242 = t206 * pkin(2) * t280 * t272;
t179 = 0.1e1 / t182;
t271 = (pkin(2) + t202) * (pkin(2) - t202);
t171 = t262 + t271;
t270 = (pkin(2) - t282) * (pkin(2) + t282);
t172 = t262 + t270;
t273 = t171 * t172;
t217 = sqrt(-t273);
t162 = 0.1e1 / t217;
t205 = 0.2e1 * t213;
t211 = t213 ^ 2;
t214 = t216 ^ 2;
t165 = t214 + (t205 - 0.2e1 / 0.3e1 * t209 - 0.2e1 / 0.3e1 * t207) * t216 + t211 + (-0.4e1 / 0.3e1 * t209 - 0.4e1 / 0.3e1 * t207) * t213 + t269 / 0.3e1;
t284 = 0.3e1 * t165;
t283 = -2 * qJD(1);
t281 = t179 / 0.2e1;
t279 = pkin(1) * t201;
t278 = pkin(2) * t200;
t258 = pkin(1) * t278;
t159 = (-t171 - t172) * t258;
t155 = qJD(1) * t159;
t276 = t155 * t162;
t275 = t159 * t162;
t186 = -pkin(2) + t279;
t274 = t162 * t186;
t268 = t200 * t217;
t267 = t201 * t216;
t266 = t202 * t282;
t265 = t216 * t200 ^ 2;
t177 = -g(1) * t278 + g(2) * t277;
t173 = g(2) * pkin(1) - t177;
t243 = g(1) * t201 + g(2) * t200;
t178 = t243 * pkin(2);
t264 = -t173 * t275 - t178 * t217;
t174 = -g(1) * pkin(1) + t178;
t263 = -t174 * t275 - t177 * t217;
t261 = -t207 / 0.2e1 + t213;
t260 = t209 - t207;
t256 = pkin(2) * t265;
t210 = 0.1e1 / pkin(3);
t255 = t210 * t272;
t195 = t201 ^ 2;
t218 = pkin(1) * t216;
t254 = (-t209 / 0.2e1 + t261) * t195 * t218;
t253 = t186 * t266;
t184 = t209 / 0.2e1 + t261;
t252 = t213 + t260;
t251 = t210 * t281;
t166 = (t213 - t209 / 0.6e1 - t207 / 0.6e1) * t216 + t211 + (-0.5e1 / 0.6e1 * t209 - 0.5e1 / 0.6e1 * t207) * t213 + t269 / 0.6e1;
t249 = -0.12e2 * t166 * t267;
t248 = 0.12e2 * t254;
t247 = t206 * t255;
t246 = t200 * t253;
t245 = qJD(1) * t255;
t244 = t180 * t210 * t258;
t241 = t179 * t242;
t240 = t162 * t247;
t239 = t162 * t245;
t238 = t162 * qJDD(1) * t255;
t237 = 0.2e1 * t241;
t236 = 0.1e1 / t172 ^ 2 * t242;
t235 = t186 * t279 - t265;
t163 = t162 / t273;
t234 = t159 * t163 * t247;
t233 = t155 * t163 * t245;
t169 = 0.1e1 / t172;
t232 = 0.1e1 / t171 ^ 2 * t169 * t242;
t231 = t162 * t210 * t241;
t230 = t235 * t217 * t266;
t175 = t252 + t262;
t229 = -0.2e1 * t256 + (-t175 * t201 - t268) * pkin(1);
t176 = t182 - t260;
t228 = 0.2e1 * t256 + (t176 * t201 - t268) * pkin(1);
t183 = t216 + t184;
t227 = 0.4e1 * t189 * t186 * t256 + (0.4e1 * pkin(1) * pkin(2) * t183 - 0.8e1 * t184 * t267) * t268 + t235 * t285;
t167 = 0.1e1 / t171;
t160 = t186 * t217;
t158 = -0.4e1 * t183 * t257 + t214 + t213 * t252 + (0.4e1 * t184 * t195 + t205 - t260) * t216;
t154 = t174 * t217;
t153 = t173 * t217;
t152 = t176 * t280 + t160;
t151 = -t175 * t280 + t160;
t150 = t152 ^ 2;
t149 = t151 ^ 2;
t148 = t159 * t274;
t147 = t155 * t274;
t1 = [qJDD(1) * MDP(1) + (g(1) * t200 - g(2) * t201) * MDP(2) + t243 * MDP(3) + (t150 * t232 + (t150 * t236 + (t150 * t237 + (-qJDD(1) * t150 + ((qJD(1) * t228 + t147) * t283 + t206 * (t148 + t228)) * t152) * t272) * t169) * t167) * MDP(4) + (-(qJD(1) * t227 + t158 * t276) * t239 + (t158 * t275 + t227) * t240 / 0.2e1 + (0.2e1 * t174 * t258 + t175 * t177 + t264) * t251 - (t174 * t175 - t153) * t244 + (-t233 - t238 + 0.2e1 * t231 + t234 / 0.2e1) * (t186 * t280 * t285 + t158 * t217)) * MDP(5) + (-0.2e1 * (-pkin(1) * t246 * t276 + (-t230 + (t249 + (pkin(1) * t284 + t248) * pkin(2)) * t200) * qJD(1)) * t239 + (-t230 + (pkin(2) * t248 + t249 + (pkin(2) * t284 - t253 * t275) * pkin(1)) * t200) * t240 + (-0.2e1 * t173 * t258 - t175 * t178 + t263) * t251 - (-t173 * t175 - t154) * t244 + (-0.2e1 * t233 - 0.2e1 * t238 + 0.4e1 * t231 + t234) * (-0.4e1 * t254 * t277 + t218 ^ 2 / 0.2e1 + (0.6e1 * t166 * t195 + 0.3e1 / 0.2e1 * t211 - t269 / 0.2e1) * t216 + (0.3e1 / 0.2e1 * t214 - t189 * t216 + t270 * t271 / 0.2e1) * t213 + (-0.3e1 * t165 * t277 - t217 * t246) * pkin(1))) * MDP(6) + (t149 * t232 + (t149 * t236 + (t149 * t237 + (-t149 * qJDD(1) + ((qJD(1) * t229 + t147) * t283 + t206 * (t148 + t229)) * t151) * t272) * t169) * t167) * MDP(7) + (((-t176 * t177 + t264) * t281 + (-t174 * t179 - (-t174 * t176 - t153) * t180) * t258) * MDP(8) + ((t176 * t178 + t263) * t281 + (t173 * t179 - (t173 * t176 - t154) * t180) * t258) * MDP(9)) / pkin(4);];
tau = t1;
