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
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1DE1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_invdynJ_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_invdynJ_fixb_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbar1DE1_invdynJ_fixb_slag_vp2: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_invdynJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE1_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1DE1_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1DE1_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:57:11
% EndTime: 2020-04-24 19:57:16
% DurationCPUTime: 2.58s
% Computational Cost: add. (4520->136), mult. (6661->268), div. (158->11), fcn. (1620->4), ass. (0->128)
t64 = sin(qJ(1));
t157 = pkin(2) * t64;
t132 = pkin(1) * t157;
t65 = cos(qJ(1));
t156 = pkin(2) * t65;
t73 = pkin(1) ^ 2;
t135 = -0.2e1 * pkin(1) * t156 + t73;
t162 = -pkin(3) - pkin(4);
t36 = (pkin(2) - t162) * (pkin(2) + t162) + t135;
t161 = -pkin(3) + pkin(4);
t37 = (pkin(2) - t161) * (pkin(2) + t161) + t135;
t141 = t36 * t37;
t74 = sqrt(-t141);
t155 = pkin(2) * t74;
t57 = pkin(1) * t65 - pkin(2);
t26 = t57 * t155;
t175 = pkin(3) ^ 2;
t134 = -pkin(4) ^ 2 + t175;
t72 = pkin(2) ^ 2;
t53 = t72 + t135;
t44 = t53 + t134;
t20 = -t132 * t44 + t26;
t178 = -t20 / 0.2e1;
t50 = 0.1e1 / t53;
t163 = t50 / 0.2e1;
t45 = t53 - t134;
t21 = t132 * t45 + t26;
t149 = Ifges(3,3) * t21;
t177 = -t149 / 0.2e1;
t121 = t64 ^ 2 * t72 * t73;
t66 = qJD(1) ^ 2;
t112 = t66 * t121;
t29 = 0.1e1 / t74;
t131 = pkin(2) * t29 * t57;
t111 = qJD(1) * t131;
t24 = (-t36 - t37) * t132;
t23 = qJD(1) * t24;
t159 = pkin(1) * t64;
t90 = qJDD(1) * t57 - t159 * t66;
t137 = t23 * t111 + t90 * t155;
t42 = (qJDD(1) * t64 + t65 * t66) * pkin(1) * pkin(2);
t176 = -t42 * t44 - 0.2e1 * t112 + t137;
t70 = 0.1e1 / pkin(3);
t117 = t70 * t163;
t119 = qJD(1) * t159;
t114 = pkin(2) * t119;
t133 = pkin(2) * qJD(1);
t118 = t74 * t133;
t25 = t57 * t118;
t18 = t114 * t45 + t25;
t51 = 0.1e1 / t53 ^ 2;
t52 = t50 * t51;
t115 = t52 * t132;
t94 = t29 * t70 * t115;
t89 = qJD(1) * t94;
t172 = t18 * t94 - 0.2e1 * t21 * t89;
t33 = 0.1e1 / t36 ^ 2;
t35 = 0.1e1 / t37 ^ 2;
t122 = t33 * t35 * t51;
t32 = 0.1e1 / t36;
t34 = 0.1e1 / t37;
t142 = t32 * t34;
t93 = t115 * t142;
t171 = -t93 + t24 * t122 / 0.2e1;
t13 = t24 * t111;
t110 = qJD(1) * t121;
t113 = pkin(1) * t118;
t126 = pkin(1) * t133;
t87 = t45 * t65 * t126 - t113 * t64 + 0.2e1 * t110;
t10 = t13 + t87;
t144 = t29 * t50;
t30 = t29 / t141;
t143 = t30 * t50;
t95 = 0.2e1 * t29 * t51 * t132;
t82 = -t143 * t24 + t95;
t169 = -t10 * t144 + t18 * t82;
t2 = t42 * t45 + 0.2e1 * t112 + t137;
t81 = qJD(1) * t95 - t143 * t23;
t168 = -t2 * t144 + t18 * t81;
t138 = t64 * t74;
t19 = (pkin(1) * t138 + t44 * t57) * pkin(2);
t16 = qJD(1) * t19;
t166 = t16 / 0.2e1;
t17 = -t114 * t44 + t25;
t165 = t17 / 0.2e1;
t164 = t19 / 0.2e1;
t158 = pkin(1) * t70;
t153 = mrSges(3,2) * t64;
t152 = mrSges(4,2) * t64;
t150 = Ifges(3,3) * t18;
t147 = t17 * t20;
t146 = t29 * t23;
t145 = t29 * t24;
t139 = t51 * t70;
t106 = 0.2e1 * t57 * t72 * t159;
t136 = qJD(1) * t106 + t65 * t113;
t46 = mrSges(4,1) * t157 + mrSges(4,2) * t156;
t48 = mrSges(3,1) * t157 + mrSges(3,2) * t156;
t124 = t29 * t139;
t123 = t51 * t142;
t116 = pkin(1) - t156;
t109 = -t124 / 0.2e1;
t108 = t124 / 0.2e1;
t105 = t123 * t150;
t102 = t21 * t23 * t30 * t139;
t101 = mrSges(3,1) * t109;
t100 = mrSges(3,2) * t108;
t99 = t18 * t109;
t98 = t18 * t108;
t84 = t18 * t89;
t80 = -0.2e1 * t110 + (-t44 * t65 - t138) * t126;
t79 = t50 * t169;
t78 = t50 * t168;
t68 = 0.1e1 / pkin(4);
t49 = (mrSges(3,1) * t65 - t153) * pkin(2);
t47 = (mrSges(4,1) * t65 - t152) * pkin(2);
t41 = mrSges(3,1) * t116 + pkin(2) * t153;
t40 = -mrSges(3,2) * pkin(1) + t48;
t39 = -mrSges(4,2) * pkin(1) + t46;
t38 = mrSges(4,1) * t116 + pkin(2) * t152;
t14 = t17 ^ 2;
t12 = t23 * t131;
t9 = t13 + t80;
t8 = t12 + t80;
t6 = (-t44 + t145) * t114 + t136;
t5 = (-qJD(1) * t44 + t146) * t132 + t136;
t3 = t66 * t106 + t42 * t74 + (t119 * t146 + t44 * t90) * pkin(2);
t1 = [t10 * t105 / 0.2e1 + Ifges(2,3) * qJDD(1) + (t78 * t177 + t79 * t150 / 0.2e1) * t29 + (t100 * t176 + t3 * t101) * t21 + t171 * Ifges(3,3) * t18 ^ 2 + t2 * t123 * t177 + (-(t145 * t41 - t44 * t49 + t48 * t74) * t117 - (m(3) * t65 + (-t40 * t50 - (-t40 * t44 + t41 * t74) * t51) * t64 * t158) * pkin(2) + mrSges(2,2) * t64 - mrSges(2,1) * t65 - ((t145 * t38 + t45 * t47 + t46 * t74) * t163 + (t39 * t50 - (t38 * t74 + t39 * t45) * t51) * t132) * t68) * g(2) + (-(t145 * t40 + t44 * t48 + t49 * t74) * t117 - (-m(3) + (t41 * t50 - (t40 * t74 + t41 * t44) * t51) * t158) * t157 + mrSges(2,1) * t64 + mrSges(2,2) * t65 - ((t145 * t39 - t45 * t46 + t47 * t74) * t163 + (-t38 * t50 - (-t38 * t45 + t39 * t74) * t51) * t132) * t68) * g(1) + (t100 * t17 + t101 * t16 - t105) * (t12 + t87) + (t102 * t165 + t8 * t98 + t9 * t99 + (-t117 * t168 - t84) * t20 + (t117 * t169 + t172) * t17) * mrSges(3,2) + (-(t165 * t9 + t166 * t6) * t51 / 0.2e1 - (-t16 ^ 2 - t14) * t115 / 0.2e1 + (-t16 * t19 - t147) * t52 * t114 + (t176 * t20 / 0.2e1 + t8 * t165 + t3 * t164 + t5 * t166) * t51 / 0.2e1) * m(3) / t175 + (t171 * t14 + ((-t144 * t9 + t17 * t82) * t165 + (-t144 * t176 + t17 * t81) * t178) * t144 + (t176 * t178 + (-t8 + t9 / 0.2e1) * t17) * t123) * Ifges(4,3) + ((t78 * t164 - t16 * t79 / 0.2e1) * t70 + t19 * t84 + t5 * t99 + t6 * t98 + (-t102 / 0.2e1 - t172) * t16) * mrSges(3,1) + (0.3e1 * qJD(1) * t93 - t23 * t122 / 0.2e1 + (t35 * t32 + t34 * t33) * t51 * t114) * (Ifges(4,3) * t147 + t18 * t149);];
tau = t1;
