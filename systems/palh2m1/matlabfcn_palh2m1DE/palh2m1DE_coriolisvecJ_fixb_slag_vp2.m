% Calculate vector of centrifugal and Coriolis load on the joints for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m1DE_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1DE_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:06
% EndTime: 2020-05-02 23:52:12
% DurationCPUTime: 1.26s
% Computational Cost: add. (1281->193), mult. (1804->286), div. (0->0), fcn. (854->16), ass. (0->132)
t101 = qJD(1) ^ 2;
t94 = sin(qJ(3));
t95 = sin(qJ(2));
t145 = t94 * t95;
t97 = cos(qJ(3));
t77 = pkin(3) * t97 + pkin(2);
t98 = cos(qJ(2));
t176 = -m(5) * (-pkin(3) * t145 + t77 * t98 + pkin(1)) - mrSges(5,1);
t183 = t101 * t176;
t133 = -qJ(4) + qJ(2);
t123 = qJ(3) + t133;
t117 = sin(t123);
t90 = qJD(2) + qJD(3);
t120 = qJD(1) + t90;
t103 = t120 * t117;
t102 = -t103 / 0.2e1;
t121 = qJD(1) - t90;
t92 = qJ(4) + qJ(2);
t86 = qJ(3) + t92;
t74 = sin(t86);
t108 = -t121 * t74 / 0.2e1;
t88 = qJD(4) + qJD(2);
t79 = qJD(3) + t88;
t159 = pkin(3) * t79;
t126 = t159 / 0.2e1;
t182 = t74 * t126 - (t102 + t108) * pkin(3);
t76 = cos(t123);
t115 = t76 * t120;
t109 = -t115 / 0.2e1;
t127 = -t159 / 0.2e1;
t75 = cos(t86);
t181 = pkin(3) * t109 - t75 * t127;
t122 = sin(t133);
t131 = qJD(1) + qJD(2);
t107 = t131 * t122;
t104 = -t107 / 0.2e1;
t132 = qJD(1) - qJD(2);
t81 = sin(t92);
t112 = -t132 * t81 / 0.2e1;
t164 = pkin(2) * t88;
t128 = t164 / 0.2e1;
t68 = -pkin(3) * t117 / 0.2e1;
t89 = -qJD(4) + qJD(2);
t80 = qJD(3) + t89;
t53 = t80 * t68;
t71 = -pkin(2) * t122 / 0.2e1;
t141 = t89 * t71 + t53;
t180 = t81 * t128 - (t104 + t112) * pkin(2) + t141 + t182;
t84 = cos(t133);
t119 = t84 * t131;
t113 = -t119 / 0.2e1;
t129 = -t164 / 0.2e1;
t83 = cos(t92);
t114 = t83 * t132 / 0.2e1;
t110 = t75 * t121 / 0.2e1;
t58 = pkin(3) * t110;
t139 = pkin(2) * t114 + t58;
t165 = pkin(2) * t84;
t160 = pkin(3) * t76;
t54 = -t80 * t160 / 0.2e1;
t140 = t54 - t89 * t165 / 0.2e1;
t179 = pkin(2) * t113 - t83 * t129 + t139 - t140 + t181;
t178 = t53 + t182;
t177 = -t54 + t58 + t181;
t175 = -0.2e1 * pkin(1);
t60 = -t97 * t98 + t145;
t48 = t60 * qJD(1);
t173 = t48 / 0.2e1;
t144 = t94 * t98;
t116 = t95 * t97 + t144;
t49 = t116 * qJD(1);
t172 = -t49 / 0.2e1;
t171 = t95 / 0.2e1;
t170 = -t98 / 0.2e1;
t163 = pkin(2) * t98;
t78 = pkin(1) + t163;
t169 = m(4) * t78;
t167 = pkin(2) * t81;
t166 = pkin(2) * t83;
t162 = pkin(3) * t74;
t161 = pkin(3) * t75;
t93 = qJ(2) + qJ(3);
t82 = sin(t93);
t158 = pkin(3) * t82;
t85 = cos(t93);
t157 = pkin(3) * t85;
t156 = pkin(3) * t90 ^ 2;
t155 = mrSges(4,3) * t48;
t154 = mrSges(4,3) * t49;
t153 = Ifges(3,4) * t95;
t152 = Ifges(3,4) * t98;
t151 = Ifges(4,4) * t49;
t150 = pkin(2) * qJD(2) ^ 2;
t148 = t77 * t95;
t147 = t78 * (-mrSges(4,1) * t49 + mrSges(4,2) * t48);
t91 = qJD(1) + qJD(4);
t146 = t91 * mrSges(6,1);
t143 = t95 * (-mrSges(4,1) * t48 - mrSges(4,2) * t49);
t99 = pkin(1) + pkin(4);
t142 = cos(qJ(4)) * t99;
t138 = t68 + t71;
t70 = t160 / 0.2e1;
t137 = t70 + t165 / 0.2e1;
t136 = pkin(2) * qJD(2);
t43 = t95 * t150 + t82 * t156;
t130 = t43 * t157;
t125 = qJD(1) * t142;
t124 = sin(qJ(4)) * t99;
t118 = qJD(1) * t124;
t111 = (Ifges(3,5) * qJD(2) + (-Ifges(3,1) * t95 - t152) * qJD(1)) * t170 + (Ifges(3,6) * qJD(2) + (-Ifges(3,2) * t98 - t153) * qJD(1)) * t171;
t17 = Ifges(4,2) * t48 + Ifges(4,6) * t90 - t151;
t40 = Ifges(4,4) * t48;
t18 = -Ifges(4,1) * t49 + Ifges(4,5) * t90 + t40;
t29 = t90 * t60;
t22 = t29 * qJD(1);
t30 = t90 * t116;
t23 = t30 * qJD(1);
t106 = t17 * t172 + t49 * (Ifges(4,1) * t48 + t151) / 0.2e1 - t90 * (Ifges(4,5) * t48 + Ifges(4,6) * t49) / 0.2e1 + Ifges(4,6) * t23 + Ifges(4,5) * t22 + t97 * t136 * t155 - (Ifges(4,2) * t49 + t18 + t40) * t48 / 0.2e1;
t105 = (Ifges(3,5) * t170 + Ifges(3,6) * t171) * qJD(2);
t69 = -t161 / 0.2e1;
t67 = -t162 / 0.2e1;
t42 = -t98 * t150 - t85 * t156;
t37 = mrSges(4,1) * t90 + t154;
t36 = -mrSges(4,2) * t90 + t155;
t9 = pkin(3) * t115 / 0.2e1 + pkin(2) * t119 / 0.2e1 + t125 + t139;
t8 = -t118 + (t103 / 0.2e1 + t108) * pkin(3) + (t107 / 0.2e1 + t112) * pkin(2);
t7 = qJD(4) * t142 + t75 * t126 + t83 * t128 + t140;
t6 = -qJD(4) * t124 + t74 * t127 + t81 * t129 + t141;
t3 = -qJD(4) * t118 + (t80 * t102 + t79 * t108) * pkin(3) + (t89 * t104 + t88 * t112) * pkin(2);
t2 = qJD(4) * t125 + (t80 * t109 + t79 * t110) * pkin(3) + (t89 * t113 + t88 * t114) * pkin(2);
t1 = t3 * mrSges(6,1);
t4 = [t90 * (Ifges(4,5) * t29 + Ifges(4,6) * t30) / 0.2e1 + t78 * (-t23 * mrSges(4,1) + t22 * mrSges(4,2)) + m(6) * (t2 * (t162 / 0.2e1 + t167 / 0.2e1 + t124 + t138) - t8 * t7 + t3 * (t161 / 0.2e1 + t166 / 0.2e1 + t142 + t137) + t9 * t6) + t6 * t146 - t42 * mrSges(5,3) + t29 * t18 / 0.2e1 + t30 * t17 / 0.2e1 + t1 + (t30 * t173 + t23 * t60) * Ifges(4,2) + (-t116 * t22 + t29 * t172) * Ifges(4,1) + (-t7 * t91 - t2) * mrSges(6,2) + (t105 + (-t143 + (-t29 * t97 + t30 * t94 + (-t116 * t94 + t60 * t97) * qJD(3)) * mrSges(4,3)) * pkin(2) + t111) * qJD(2) + (-t116 * t23 + t30 * t172 + t29 * t173 + t22 * t60) * Ifges(4,4) + (t78 * (-mrSges(4,1) * t30 + mrSges(4,2) * t29) - 0.2e1 * t176 * (-qJD(2) * t148 + (-qJD(2) * t144 - t116 * qJD(3)) * pkin(3)) + ((mrSges(3,2) * t175 + 0.3e1 / 0.2e1 * t152) * t98 + (mrSges(3,1) * t175 - 0.3e1 / 0.2e1 * t153 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t98 + (mrSges(4,1) * t60 + mrSges(4,2) * t116 - 0.2e1 * t169) * pkin(2)) * t95) * qJD(2)) * qJD(1); (mrSges(6,1) * t180 + mrSges(6,2) * t179) * t91 + (t101 * t95 * t169 + (-t22 * t97 + (-qJD(2) * t49 + t23) * t94) * mrSges(4,3) + (t36 * t97 - t37 * t94 + (-mrSges(4,1) * t94 - mrSges(4,2) * t97) * qJD(2)) * qJD(3)) * pkin(2) + t42 * (-pkin(2) * t95 - t158) * m(5) + (-t147 + pkin(2) * t143 + ((-Ifges(3,1) * t98 + t153) * t171 + t98 * (Ifges(3,2) * t95 - t152) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t95 - mrSges(3,2) * t98)) * qJD(1) + t105 - t111) * qJD(1) + t106 + (-pkin(3) * t144 - t148) * t183 + (t2 * (t67 - t167 / 0.2e1 + t138) + t3 * (t69 - t166 / 0.2e1 + t137) + t180 * t9 + t179 * t8) * m(6) + (-m(6) - m(5)) * t43 * (t157 + t163); ((-mrSges(4,2) * qJD(3) - t36) * t97 + (-mrSges(4,1) * qJD(3) - t154 + t37) * t94) * t136 + (mrSges(6,1) * t178 + mrSges(6,2) * t177) * t91 + m(5) * (-t42 * t158 - t130) + t106 - qJD(1) * t147 - t116 * pkin(3) * t183 + (t2 * (t68 + t67) + t3 * (t70 + t69) - t130 + t178 * t9 + t177 * t8) * m(6); -t8 * t146 + t1 + (t9 * t91 - t2) * mrSges(6,2);];
tauc = t4(:);
