% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% fivebar1OL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = fivebar1OL_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:17
% EndTime: 2020-04-27 06:13:18
% DurationCPUTime: 0.77s
% Computational Cost: add. (296->126), mult. (456->135), div. (0->0), fcn. (276->8), ass. (0->80)
t153 = qJD(3) + qJD(4);
t149 = t153 ^ 2;
t135 = t149 * pkin(1) - g(1);
t158 = sin(qJ(3));
t162 = cos(qJ(3));
t151 = -qJDD(4) - qJDD(3);
t175 = t151 * pkin(1) + g(2);
t197 = -t135 * t162 + t175 * t158;
t157 = sin(qJ(4));
t196 = -0.2e1 * t157;
t159 = sin(qJ(2));
t195 = 0.2e1 * t159;
t194 = pkin(1) * g(2);
t193 = pkin(2) * (qJDD(1) + qJDD(2) / 0.2e1);
t192 = pkin(3) * (qJDD(3) + qJDD(4) / 0.2e1);
t161 = cos(qJ(4));
t191 = (t157 * t162 + t158 * t161) * g(3);
t190 = (-t157 * t158 + t161 * t162) * g(3);
t160 = sin(qJ(1));
t163 = cos(qJ(2));
t164 = cos(qJ(1));
t189 = (-t160 * t159 + t163 * t164) * g(3);
t188 = t157 * g(3);
t187 = t158 * g(1);
t186 = t158 * g(3);
t185 = t160 * g(1);
t184 = t160 * g(3);
t183 = t161 * g(3);
t182 = t162 * g(3);
t181 = t163 * g(3);
t180 = t164 * g(3);
t179 = qJDD(3) * pkin(3);
t142 = t158 * g(2);
t146 = t162 * g(1);
t130 = t146 + t142;
t144 = t160 * g(2);
t148 = t164 * g(1);
t131 = t148 + t144;
t177 = qJD(2) * pkin(2) * (qJD(1) + qJD(2) / 0.2e1);
t176 = qJD(4) * pkin(3) * (qJD(3) + qJD(4) / 0.2e1);
t145 = t162 * g(2);
t128 = -t145 + t187;
t147 = t164 * g(2);
t129 = -t147 + t185;
t154 = qJD(1) + qJD(2);
t150 = t154 ^ 2;
t152 = qJDD(1) + qJDD(2);
t174 = t159 * t150 - t152 * t163;
t110 = t157 * t149 + t151 * t161;
t113 = t149 * t161 - t157 * t151;
t173 = t110 * t158 - t113 * t162;
t114 = t150 * t163 + t159 * t152;
t172 = -t114 * t164 + t160 * t174;
t120 = t128 + t179;
t165 = qJD(3) ^ 2;
t122 = -t165 * pkin(3) - t130;
t171 = t161 * t120 - t157 * t122;
t170 = t157 * t120 + t161 * t122;
t121 = qJDD(1) * pkin(2) + t129;
t166 = qJD(1) ^ 2;
t123 = -t166 * pkin(2) - t131;
t169 = t163 * t121 - t159 * t123;
t168 = t159 * t121 + t163 * t123;
t136 = pkin(1) * qJDD(3) - g(2);
t137 = pkin(1) * t165 - g(1);
t167 = t136 * t162 - t158 * t137;
t127 = t164 * qJDD(1) - t160 * t166;
t125 = -t160 * qJDD(1) - t164 * t166;
t126 = t162 * qJDD(3) - t158 * t165;
t124 = -t158 * qJDD(3) - t162 * t165;
t143 = t159 * g(3);
t132 = -0.2e1 * t176;
t116 = pkin(2) * t121;
t115 = (t159 * t164 + t160 * t163) * g(3);
t108 = -t175 * t162 - t135 * t158 + pkin(3) * (qJDD(4) + 0.2e1 * qJDD(3));
t107 = t160 * t114 + t164 * t174;
t106 = -t110 * t162 - t158 * t113;
t105 = (0.2e1 * t177 - t131) * t163 + (-t147 / 0.2e1 + t185 / 0.2e1 + t193) * t195;
t104 = (-t129 - 0.2e1 * t193) * t163 + (-t148 / 0.2e1 - t144 / 0.2e1 + t177) * t195;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t127, 0, t125, 0, -t184, -t180, g(2), 0, 0, 0, t107, 0, -t172, 0, t115, t189, -t127 * pkin(2) + g(2), -pkin(2) * t184, 0, 0, t126, 0, t124, 0, -t186, -t182, g(2), 0, 0, 0, t106, 0, t173, 0, -t191, -t190, -t126 * pkin(3) + g(2), -pkin(3) * t186; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t125, 0, t127, 0, t180, -t184, -g(1), 0, 0, 0, t172, 0, t107, 0, -t189, t115, t125 * pkin(2) - g(1), pkin(2) * t180, 0, 0, -t124, 0, t126, 0, t182, -t186, -g(1), pkin(1) * g(3), 0, 0, -t173, 0, t106, 0, t190, -t191, t124 * pkin(3) - g(1), (pkin(3) * t162 + pkin(1)) * g(3); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t129, t131, 0, 0, 0, 0, 0, 0, 0, t152, t104, t105, 0, t116, 0, 0, 0, 0, 0, qJDD(3), t167, -t158 * t136 - t137 * t162, 0, -t194, 0, 0, 0, 0, 0, -t151, t108 * t161 - (0.2e1 * t176 - t197) * t157, (t132 + t197) * t161 - t157 * t108, 0, -t194 + (t167 + t179) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t166, 0, 0, -g(3), -t129, 0, 0, 0, t174, 0, t114, 0, t143, t181, -t121, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, 0, qJDD(1), 0, g(3), 0, -t131, 0, 0, 0, -t114, 0, t174, 0, -t181, t143, t123, pkin(2) * g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t129, t131, 0, 0, 0, 0, 0, 0, 0, t152, t104, t105, 0, t116, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, 0, -t150, 0, 0, -g(3), t169, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, 0, t152, 0, g(3), 0, -t168, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, -t169, t168, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t165, 0, 0, -g(3), -t128, 0, 0, 0, -t110, 0, -t113, 0, -t188, -t183, -t120, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, 0, qJDD(3), 0, g(3), 0, -t130, 0, 0, 0, t113, 0, -t110, 0, t183, -t188, t122, pkin(3) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t128, t130, 0, 0, 0, 0, 0, 0, 0, -t151, (t128 + 0.2e1 * t192) * t161 + (-t146 / 0.2e1 - t142 / 0.2e1 + t176) * t196, (t132 + t130) * t161 + (-t145 / 0.2e1 + t187 / 0.2e1 + t192) * t196, 0, pkin(3) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, 0, -t149, 0, 0, -g(3), -t171, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, 0, -t151, 0, g(3), 0, t170, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t171, -t170, 0, 0;];
m_new_reg = t1;
