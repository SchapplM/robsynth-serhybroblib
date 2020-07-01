% Calculate vector of inverse dynamics joint torques for
% fourbarprisDE2
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisDE2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:23
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbarprisDE2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_invdynJ_fixb_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_invdynJ_fixb_mdp_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbarprisDE2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisDE2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_invdynJ_fixb_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbarprisDE2_invdynJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:23:20
% EndTime: 2020-06-27 17:23:23
% DurationCPUTime: 1.04s
% Computational Cost: add. (1252->79), mult. (685->158), div. (184->15), fcn. (38->2), ass. (0->78)
t143 = (qJ(1) + pkin(3));
t183 = -pkin(2) - t143;
t135 = pkin(1) + t183;
t198 = 1 / t135;
t208 = 2 * t198;
t132 = pkin(1) - t183;
t182 = -pkin(2) + t143;
t133 = pkin(1) - t182;
t134 = pkin(1) + t182;
t193 = t134 * t135;
t177 = t133 * t193;
t169 = t132 * t177;
t152 = sqrt(-t169);
t114 = 0.1e1 / t152;
t207 = (-t177 + (t193 + (t134 - t135) * t133) * t132) * t114;
t206 = 2 * g(1);
t205 = 0.1e1 / t169 * t207;
t122 = 0.1e1 / t132;
t126 = 0.1e1 / t134;
t204 = t122 * t126;
t139 = (t143 ^ 2);
t140 = 1 / t143;
t146 = (qJ(1) ^ 2);
t147 = (pkin(3) ^ 2);
t195 = (pkin(3) * qJ(1));
t176 = -2 * t195 - t146 - t147;
t148 = pkin(2) ^ 2;
t150 = (pkin(1) ^ 2);
t186 = (t150 - t148);
t119 = t176 + t186;
t200 = t119 ^ 2;
t203 = t140 / t139 * t200;
t199 = 0.1e1 / t133;
t202 = t198 * t199;
t184 = qJDD(1) * t200;
t201 = t200 * t202;
t125 = 0.1e1 / t133 ^ 2;
t129 = 0.1e1 / t135 ^ 2;
t151 = 0.1e1 / pkin(1);
t197 = -t151 / 0.2e1;
t196 = t151 / 0.2e1;
t145 = (qJD(1) ^ 2);
t194 = t119 * t145;
t192 = t139 * t145;
t141 = 1 / (t143 ^ 2);
t191 = t141 * t145;
t190 = t143 * t145;
t112 = t152 * g(2);
t118 = -t176 + t186;
t175 = t207 / 0.2e1;
t107 = g(2) * t175;
t161 = ((t143 * t206) + t107) * t196;
t174 = t141 * t197;
t189 = t140 * t161 + ((g(1) * t118) + t112) * t174;
t188 = (t147 * g(1)) + t112;
t187 = (-t148 - t150);
t121 = t143 * qJD(1);
t185 = qJD(1) * t121;
t181 = 4 * t185;
t180 = t200 * t191;
t179 = t145 * t203;
t178 = t198 * t192;
t173 = t198 * t180;
t172 = t129 * t180;
t170 = t126 * t178;
t168 = t126 * t173;
t127 = 0.1e1 / t134 ^ 2;
t167 = t127 * t173;
t166 = t146 * t202 * t204;
t165 = t199 * t168;
t164 = t125 * t168;
t163 = t141 * t166;
t162 = t145 * t166;
t106 = g(1) * t175 - 0.2e1 * g(2) * t143;
t123 = 0.1e1 / t132 ^ 2;
t160 = t123 * t165;
t111 = g(1) * t152 - g(2) * t118;
t1 = [(t160 / 0.2e1 + (-t164 / 0.2e1 + (t167 / 0.2e1 + (-t172 / 0.2e1 + ((t179 + (-t184 + 2 * (2 * t185 - t190) * t119) * t141) * t198)) * t126) * t199) * t122) * MDP(1) + t189 * MDP(2) + (-t106 * t140 / 0.2e1 + t111 * t141 / 0.2e1) * t151 * MDP(3) + ((t111 * t197 - t114 * t194) * t141 + (-t194 * t205 / 0.2e1 + t106 * t196 + 0.2e1 * (qJDD(1) * t119 + t190) * t114 + (qJD(1) * t119 * t205 - 0.4e1 * t114 * t121) * qJD(1)) * t140) * MDP(4) + (-t122 * t165 + (t160 + (-t164 + (t167 + (-t172 + ((t179 + (-t184 + (t181 - 2 * t190) * t119) * t141) * t208)) * t126) * t199) * t122) * qJ(1) + t189) * MDP(5) + (-t163 * t184 + t162 * t203 + qJDD(1) + (qJ(1) * t107 + ((3 * t146 + t187 + 4 * t195) * g(1)) + t188) * t140 * t196 + (((t146 - t150) * pkin(3) * t206) + (t188 + ((t146 + t187) * g(1))) * qJ(1)) * t174 + (-0.2e1 * t140 * t162 + t163 * t181) * t119 + ((t122 * t127 + t123 * t126) * t146 * t201 / 0.2e1 + (-qJ(1) * t201 - (t125 * t198 + t129 * t199) * t200 * t146 / 0.2e1) * t204) * t191) * MDP(6) + 0.2e1 * (t123 * t199 * t170 + (-t125 * t170 + (t127 * t178 + (-t129 * t192 + ((-qJDD(1) * t139 - t190) * t208)) * t126) * t199) * t122) * MDP(7) + (MDP(9) * t106 * t197 + MDP(8) * t161) / pkin(2);];
tau = t1;
