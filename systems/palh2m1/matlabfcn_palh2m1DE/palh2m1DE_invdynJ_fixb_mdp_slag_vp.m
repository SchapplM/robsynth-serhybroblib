% Calculate vector of inverse dynamics joint torques for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m1DE_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m1DE_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'palh2m1DE_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:37
% EndTime: 2020-06-30 17:39:39
% DurationCPUTime: 1.09s
% Computational Cost: add. (525->153), mult. (1092->238), div. (0->0), fcn. (803->8), ass. (0->90)
t180 = cos(qJ(2));
t161 = pkin(2) * t180 + pkin(1);
t206 = 2 * qJD(1);
t243 = t161 * t206;
t177 = sin(qJ(1));
t181 = cos(qJ(1));
t152 = g(1) * t177 - g(2) * t181;
t176 = sin(qJ(2));
t218 = qJD(2) * t176;
t237 = 2 * qJDD(1);
t242 = pkin(2) * t206 * t218 - t161 * t237 - t152;
t168 = qJD(2) + qJD(3);
t175 = sin(qJ(3));
t225 = t180 * t175;
t179 = cos(qJ(3));
t226 = t179 * t176;
t190 = t225 + t226;
t134 = t168 * t190;
t160 = pkin(3) * t179 + pkin(2);
t229 = t176 * t175;
t147 = t179 * t180 - t229;
t188 = t147 * qJD(3);
t217 = qJD(2) * t180;
t129 = t160 * t217 + (-t175 * t218 + t188) * pkin(3);
t189 = t147 * qJD(2);
t133 = t188 + t189;
t142 = pkin(3) * t225 + t160 * t176;
t153 = g(1) * t181 + g(2) * t177;
t241 = qJD(2) * t129 + qJDD(2) * t142 + (qJD(3) * t133 + qJDD(3) * t190) * pkin(3) + t153;
t193 = g(3) * t175 - t153 * t179;
t200 = t179 * g(3) + t175 * t153;
t183 = qJD(1) ^ 2;
t234 = t161 * t183;
t240 = t179 * qJDD(2) * pkin(2) - t193 * t176 + t180 * t200 + t190 * t234;
t151 = -0.2e1 * t225 * t226;
t172 = t179 ^ 2;
t173 = t180 ^ 2;
t211 = t173 - 0.1e1 / 0.2e1;
t239 = t151 + (-t175 ^ 2 + t172) * t211;
t163 = t172 - 0.1e1 / 0.2e1;
t171 = t176 ^ 2;
t221 = t171 - t173;
t238 = -t221 * t163 + t151;
t235 = pkin(3) * qJD(3);
t167 = qJDD(1) + qJDD(4);
t174 = sin(qJ(4));
t233 = t167 * t174;
t178 = cos(qJ(4));
t232 = t167 * t178;
t169 = qJD(1) + qJD(4);
t231 = t169 * t174;
t230 = t169 * t178;
t228 = t176 * t180;
t227 = t176 * t183;
t224 = t241 * t178;
t143 = t147 * t235;
t195 = -pkin(3) * t229 + t160 * t180;
t223 = -qJD(2) * t195 + t129 - t143;
t220 = MDP(11) * t190;
t144 = t190 * pkin(3);
t219 = qJD(1) * t144;
t216 = qJD(4) * t178;
t128 = -t160 * t218 + (-qJD(3) * t190 - t175 * t217) * pkin(3);
t215 = t128 * qJD(1);
t214 = t183 * MDP(18);
t213 = -qJD(1) - t169;
t212 = -qJD(4) + t169;
t210 = qJD(1) * qJD(2);
t209 = qJDD(1) * t190;
t208 = qJDD(1) * t180;
t191 = pkin(1) + t195;
t139 = pkin(4) + t191;
t207 = t139 * qJDD(1);
t205 = -2 * MDP(12) * t183;
t201 = -0.2e1 * pkin(1) * t210;
t199 = t212 * t139;
t198 = qJD(2) * (-qJD(3) + t168);
t197 = qJD(3) * (-qJD(2) - t168);
t196 = qJD(4) * t213;
t194 = t147 * t234 - t176 * t200 - t180 * t193;
t131 = qJD(2) * t142 + t190 * t235;
t192 = t131 * t216 + (t207 + t215) * t178 + t241 * t174;
t187 = pkin(1) * t183 + t153;
t186 = pkin(1) * t237 + t152;
t166 = qJDD(2) + qJDD(3);
t185 = -t183 * t147 * t220 + (-t209 + (t147 * t168 - t133) * qJD(1)) * MDP(13) - qJDD(1) * t147 * MDP(14) + t166 * MDP(15);
t182 = qJD(2) ^ 2;
t159 = t167 * MDP(19);
t132 = pkin(3) * t189 + t143;
t1 = [(qJDD(1) * MDP(1)) + t152 * MDP(2) + t153 * MDP(3) + (qJDD(1) * t171 + 0.2e1 * t210 * t228) * MDP(4) + 0.2e1 * (t176 * t208 - t221 * t210) * MDP(5) + (-qJDD(2) * t176 - t180 * t182) * MDP(6) + (-qJDD(2) * t180 + t176 * t182) * MDP(7) + (t176 * t201 + t180 * t186) * MDP(9) + (-t176 * t186 + t180 * t201) * MDP(10) - (-0.2e1 * qJD(1) * t133 - t209) * t220 + 0.4e1 * ((t211 * t179 * t175 + t163 * t228) * qJDD(1) + t238 * t210 + t239 * qJD(1) * qJD(3)) * MDP(12) + (-t133 * t168 - t166 * t190) * MDP(13) + (t134 * t168 - t147 * t166) * MDP(14) + (-t134 * t243 - t242 * t147) * MDP(16) + (-t133 * t243 + t242 * t190) * MDP(17) + (t191 * t237 + t152 + 0.2e1 * t215) * MDP(18) + t159 + ((t128 * t169 + t152) * t178 + (t174 * t196 + t232) * t139 + t192) * MDP(20) + (t178 * t139 * t196 + (-qJD(4) * t131 - t152 + (-qJDD(1) - t167) * t139 + t213 * t128) * t174 + t224) * MDP(21); -t180 * MDP(4) * t227 + t221 * t183 * MDP(5) - t176 * qJDD(1) * MDP(6) - MDP(7) * t208 + qJDD(2) * MDP(8) + (g(3) * t180 + t176 * t187) * MDP(9) + (-g(3) * t176 + t180 * t187) * MDP(10) + t238 * t205 + ((t147 * t227 + t179 * t166 + t175 * t197) * pkin(2) + t240) * MDP(16) + ((-t190 * t227 + (-qJDD(2) - t166) * t175 + t179 * t197) * pkin(2) + t194) * MDP(17) + t142 * t214 + (t142 * t233 + (t142 * t230 + t223 * t174) * t169) * MDP(20) + (t142 * t232 + (-t142 * t231 + t223 * t178) * t169) * MDP(21) + t185; t239 * t205 + (pkin(2) * t175 * t198 + t240) * MDP(16) + ((-qJDD(2) * t175 + t179 * t198) * pkin(2) + t194) * MDP(17) + t144 * t214 + (-(t132 * t174 - t178 * t219) * t169 + (t133 * t231 - (-t169 * t216 - t233) * t190) * pkin(3)) * MDP(20) + ((-t132 * t178 - t174 * t219) * t169 + (t133 * t230 - (qJD(4) * t231 - t232) * t190) * pkin(3)) * MDP(21) + t185; t159 + t192 * MDP(20) + t224 * MDP(21) + ((-t131 * t169 + t152) * MDP(20) + MDP(21) * qJD(1) * t199) * t178 + ((t212 * t131 - t152 - t207) * MDP(21) + (MDP(20) * t199 - t128 * MDP(21)) * qJD(1)) * t174;];
tau = t1;
