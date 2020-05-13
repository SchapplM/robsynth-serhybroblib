% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbar1TE
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
% tauc [1x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1TE_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_coriolisvecJ_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_coriolisvecJ_fixb_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_coriolisvecJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1TE_coriolisvecJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1TE_coriolisvecJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1TE_coriolisvecJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:49:17
% EndTime: 2020-04-24 19:49:21
% DurationCPUTime: 1.40s
% Computational Cost: add. (4139->94), mult. (6052->198), div. (142->10), fcn. (1454->4), ass. (0->110)
t138 = pkin(1) * pkin(2);
t48 = sin(qJ(1));
t110 = t48 * t138;
t49 = cos(qJ(1));
t131 = pkin(1) * t49;
t45 = -pkin(2) + t131;
t130 = pkin(2) * t45;
t56 = pkin(1) ^ 2;
t113 = -0.2e1 * pkin(2) * t131 + t56;
t134 = -pkin(3) - pkin(4);
t34 = (pkin(2) - t134) * (pkin(2) + t134) + t113;
t133 = -pkin(3) + pkin(4);
t35 = (pkin(2) - t133) * (pkin(2) + t133) + t113;
t118 = t34 * t35;
t57 = sqrt(-t118);
t24 = t57 * t130;
t145 = pkin(3) ^ 2;
t112 = -pkin(4) ^ 2 + t145;
t55 = pkin(2) ^ 2;
t42 = t55 + t113;
t37 = t42 + t112;
t19 = -t110 * t37 + t24;
t146 = -t19 / 0.2e1;
t111 = pkin(2) * qJD(1);
t97 = t57 * t111;
t23 = t45 * t97;
t102 = pkin(1) * t111;
t94 = t48 * t102;
t16 = -t37 * t94 + t23;
t136 = t16 / 0.2e1;
t38 = t42 - t112;
t17 = t38 * t94 + t23;
t20 = t110 * t38 + t24;
t27 = 0.1e1 / t57;
t53 = 0.1e1 / pkin(3);
t39 = 0.1e1 / t42;
t40 = 0.1e1 / t42 ^ 2;
t41 = t39 * t40;
t96 = t41 * t110;
t76 = t27 * t53 * t96;
t73 = qJD(1) * t76;
t144 = t17 * t76 - 0.2e1 * t20 * t73;
t22 = (-t34 - t35) * t110;
t30 = 0.1e1 / t34;
t32 = 0.1e1 / t35;
t119 = t30 * t32;
t75 = t96 * t119;
t31 = 0.1e1 / t34 ^ 2;
t33 = 0.1e1 / t35 ^ 2;
t99 = t31 * t33 * t40;
t142 = -t75 + t22 * t99 / 0.2e1;
t121 = t27 * t39;
t28 = t27 / t118;
t120 = t28 * t39;
t77 = 0.2e1 * t27 * t40 * t110;
t65 = -t120 * t22 + t77;
t109 = t27 * t130;
t92 = qJD(1) * t109;
t12 = t22 * t92;
t117 = t38 * t49;
t98 = t48 ^ 2 * t55 * t56;
t91 = qJD(1) * t98;
t93 = pkin(1) * t97;
t70 = t102 * t117 - t48 * t93 + 0.2e1 * t91;
t8 = t12 + t70;
t140 = -t8 * t121 + t17 * t65;
t21 = qJD(1) * t22;
t10 = t21 * t92;
t115 = t48 * t57;
t50 = qJD(1) ^ 2;
t5 = t10 + (0.2e1 * t98 + (-t115 + t117) * t138) * t50;
t64 = qJD(1) * t77 - t120 * t21;
t139 = -t5 * t121 + t17 * t64;
t18 = (pkin(1) * t115 + t37 * t45) * pkin(2);
t15 = t18 * qJD(1);
t137 = t15 / 0.2e1;
t135 = t18 / 0.2e1;
t72 = (-t37 * t49 - t115) * t138;
t4 = t10 + (t72 - 0.2e1 * t98) * t50;
t129 = t4 * t19;
t126 = Ifges(3,3) * t17;
t125 = Ifges(3,3) * t20;
t123 = t16 * t19;
t122 = t27 * t21;
t116 = t40 * t53;
t95 = 0.2e1 * t45 * t48 * t55;
t114 = qJD(1) * pkin(1) * t95 + t49 * t93;
t101 = t27 * t116;
t100 = t40 * t119;
t90 = -t101 / 0.2e1;
t89 = t101 / 0.2e1;
t88 = -t100 / 0.2e1;
t85 = t100 * t126;
t84 = t20 * t21 * t28 * t116;
t83 = mrSges(3,1) * t90;
t82 = mrSges(3,2) * t89;
t81 = t17 * t90;
t80 = t17 * t89;
t67 = t17 * t73;
t63 = qJD(1) * t72 - 0.2e1 * t91;
t62 = t39 * t140;
t61 = t39 * t139;
t13 = t16 ^ 2;
t11 = t21 * t109;
t9 = t12 + t63;
t6 = t11 + t63;
t3 = (t27 * t22 - t37) * t94 + t114;
t2 = (-qJD(1) * t37 + t122) * t110 + t114;
t1 = (t50 * t95 + (t50 * t49 * t57 + (qJD(1) * t122 - t50 * t37) * t48) * pkin(2)) * pkin(1);
t7 = [t8 * t85 / 0.2e1 + t27 * t62 * t126 / 0.2e1 + (t1 * t83 + t4 * t82) * t20 + t142 * Ifges(3,3) * t17 ^ 2 + (t5 * t88 - t27 * t61 / 0.2e1) * t125 + (t15 * t83 + t16 * t82 - t85) * (t11 + t70) + (-(t136 * t9 + t137 * t3) * t40 / 0.2e1 - (-t15 ^ 2 - t13) * t96 / 0.2e1 + (-t15 * t18 - t123) * t41 * t94 + (t129 / 0.2e1 + t6 * t136 + t1 * t135 + t2 * t137) * t40 / 0.2e1) * m(3) / t145 + ((t140 * t136 + t139 * t146) * t39 * t53 + t136 * t84 + t144 * t16 - t19 * t67 + t6 * t80 + t81 * t9) * mrSges(3,2) + ((t61 * t135 - t15 * t62 / 0.2e1) * t53 + t18 * t67 + t2 * t81 + t3 * t80 + (-t84 / 0.2e1 - t144) * t15) * mrSges(3,1) + (t88 * t129 + (-t6 + t9 / 0.2e1) * t16 * t100 + t142 * t13 + ((-t121 * t9 + t16 * t65) * t136 + (-t121 * t4 + t16 * t64) * t146) * t121) * Ifges(4,3) + (0.3e1 * qJD(1) * t75 - t21 * t99 / 0.2e1 + (t33 * t30 + t32 * t31) * t40 * t94) * (Ifges(4,3) * t123 + t17 * t125);];
tauc = t7(:);
