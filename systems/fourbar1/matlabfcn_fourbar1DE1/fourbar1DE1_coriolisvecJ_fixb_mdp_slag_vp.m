% Calculate minimal parameter regressor of Coriolis joint torque vector for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1DE1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [1x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:30
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1DE1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbar1DE1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:29:54
% EndTime: 2020-06-26 17:30:02
% DurationCPUTime: 1.08s
% Computational Cost: add. (1342->91), mult. (2206->175), div. (50->9), fcn. (396->4), ass. (0->81)
t173 = pkin(1) ^ 2;
t160 = cos(qJ(1));
t219 = pkin(2) * t160;
t198 = pkin(1) * t219;
t202 = -0.2e1 * t198 + t173;
t170 = pkin(2) ^ 2;
t142 = t170 + t202;
t165 = pkin(4) ^ 2;
t166 = pkin(3) ^ 2;
t148 = t165 + t166;
t161 = pkin(3) + pkin(4);
t222 = pkin(4) - pkin(3);
t207 = (t222 ^ 2 * t161 ^ 2);
t229 = 0.2e1 * t142 * t148 - (2 * t207);
t159 = sin(qJ(1));
t221 = pkin(1) * t159;
t228 = pkin(2) * t221;
t220 = pkin(1) * t160;
t146 = -pkin(2) + t220;
t209 = (pkin(2) + t161) * (pkin(2) - t161);
t136 = t202 + t209;
t208 = (pkin(2) - t222) * (pkin(2) + t222);
t137 = t202 + t208;
t211 = t136 * t137;
t174 = sqrt(-t211);
t125 = t146 * t174;
t200 = t166 - t165;
t196 = t170 + t200;
t138 = t196 + t202;
t120 = -t138 * t221 + t125;
t139 = t142 - t200;
t121 = t139 * t221 + t125;
t227 = MDP(4) * t121 ^ 2 + MDP(7) * t120 ^ 2;
t127 = 0.1e1 / t174;
t226 = t127 ^ 2;
t224 = -2 * qJD(1);
t223 = pkin(1) * pkin(2);
t216 = t127 * t146;
t201 = -t165 / 0.2e1 + t170;
t144 = t166 / 0.2e1 + t201;
t143 = t173 + t144;
t154 = t160 ^ 2;
t163 = 0.2e1 * t170;
t171 = t173 ^ 2;
t123 = -0.4e1 * t143 * t198 + t171 + t170 * t196 + (0.4e1 * t144 * t154 + t163 - t200) * t173;
t215 = t226 * t123;
t168 = t170 ^ 2;
t130 = t171 + (t163 - 0.2e1 / 0.3e1 * t166 - 0.2e1 / 0.3e1 * t165) * t173 + t168 + (-0.4e1 / 0.3e1 * t166 - 0.4e1 / 0.3e1 * t165) * t170 + t207 / 0.3e1;
t131 = (t170 - t166 / 0.6e1 - t165 / 0.6e1) * t173 + t168 + (-0.5e1 / 0.6e1 * t166 - 0.5e1 / 0.6e1 * t165) * t170 + t207 / 0.6e1;
t175 = pkin(1) * t173;
t204 = t161 * t222;
t210 = t146 * t159;
t193 = t204 * t210;
t194 = pkin(2) * (-t166 / 0.2e1 + t201) * t154 * t175;
t114 = -0.4e1 * t160 * t194 + t175 ^ 2 / 0.2e1 + (0.6e1 * t131 * t154 + 0.3e1 / 0.2e1 * t168 - t207 / 0.2e1) * t173 + (0.3e1 / 0.2e1 * t171 - t148 * t173 + t208 * t209 / 0.2e1) * t170 + (-0.3e1 * t130 * t219 - t174 * t193) * pkin(1);
t128 = t127 / t211;
t214 = t128 * t114;
t115 = pkin(1) * t210 * t229 + t123 * t174;
t213 = t128 * t115;
t132 = 0.1e1 / t136;
t134 = 0.1e1 / t137;
t212 = t132 * t134;
t206 = t159 * t174;
t205 = t160 * t173;
t203 = t173 * t159 ^ 2;
t197 = pkin(2) * t203;
t124 = (-t136 - t137) * t228;
t191 = t226 * pkin(1) * t193;
t190 = t146 * t220 - t203;
t189 = (0.4e1 * t143 * t223 - 0.8e1 * t144 * t205) * t174 * MDP(5);
t188 = -0.2e1 * t197 + (-t138 * t160 - t206) * pkin(1);
t187 = 0.2e1 * t197 + (t139 * t160 - t206) * pkin(1);
t186 = t190 * MDP(6) * t174 * t204;
t185 = (0.3e1 * t130 * t223 - 0.12e2 * t131 * t205 + 0.12e2 * t194) * MDP(6);
t184 = (0.4e1 * t146 * t148 * t197 + t190 * t229) * MDP(5);
t167 = 0.1e1 / pkin(3);
t164 = qJD(1) ^ 2;
t122 = qJD(1) * t124;
t117 = t124 * t216;
t116 = t122 * t216;
t1 = [((((qJD(1) * t187 + t116) * t224 + t164 * (t117 + t187)) * MDP(4) * t121 + ((qJD(1) * t188 + t116) * t224 + t164 * (t117 + t188)) * MDP(7) * t120) * t212 + (((t184 / 0.2e1 - t186 + (t189 / 0.2e1 + t185) * t159) * t127 + 0.2e1 * ((t213 / 0.4e1 + t215 / 0.4e1) * MDP(5) + (t214 / 0.2e1 - t191 / 0.2e1) * MDP(6)) * t124) * t164 + ((-t184 + 0.2e1 * t186 + (-t189 - 0.2e1 * t185) * t159) * t127 * qJD(1) + 0.2e1 * ((-t213 / 0.2e1 - t215 / 0.2e1) * MDP(5) + (t191 - t214) * MDP(6)) * t122) * qJD(1)) * t167 + ((0.2e1 * t227 * t212 + (0.2e1 * MDP(5) * t115 + 0.4e1 * MDP(6) * t114) * t167 * t127) / t142 + t227 * (t132 / t137 ^ 2 + 0.1e1 / t136 ^ 2 * t134)) * t164 * t228) * t170 / t142 ^ 2;];
tauc = t1;
