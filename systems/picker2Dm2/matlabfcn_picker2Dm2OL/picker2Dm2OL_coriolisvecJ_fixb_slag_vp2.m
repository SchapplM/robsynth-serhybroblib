% Calculate vector of centrifugal and Coriolis load on the joints for
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [12x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = picker2Dm2OL_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_coriolisvecJ_fixb_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_coriolisvecJ_fixb_slag_vp2: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_coriolisvecJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2OL_coriolisvecJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm2OL_coriolisvecJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:55
% EndTime: 2020-05-09 23:19:03
% DurationCPUTime: 1.82s
% Computational Cost: add. (1760->216), mult. (4128->341), div. (0->0), fcn. (2176->14), ass. (0->136)
t164 = mrSges(9,2) * cos(qJ(8)) + mrSges(9,1) * sin(qJ(8));
t88 = sin(qJ(2));
t94 = cos(qJ(2));
t163 = mrSges(3,1) * t88 + mrSges(3,2) * t94;
t162 = t164 * pkin(1) * qJD(8);
t80 = qJD(1) + qJD(2);
t75 = qJD(6) + t80;
t136 = t75 * mrSges(7,1);
t143 = mrSges(7,2) * t75;
t85 = sin(qJ(6));
t91 = cos(qJ(6));
t103 = t85 * t88 - t91 * t94;
t130 = pkin(1) * qJD(1);
t52 = t103 * t130;
t104 = t85 * t94 + t88 * t91;
t55 = t104 * t130;
t161 = -t55 * t136 + t143 * t52;
t116 = t88 * t130;
t74 = t94 * t130;
t62 = pkin(2) * t80 + t74;
t87 = sin(qJ(3));
t93 = cos(qJ(3));
t49 = t116 * t93 + t62 * t87;
t83 = sin(qJ(9));
t139 = t49 * t83;
t64 = t87 * t116;
t47 = -t62 * t93 + t64;
t77 = qJD(3) + t80;
t44 = pkin(6) * t77 + t47;
t89 = cos(qJ(9));
t14 = -t44 * t89 - t139;
t125 = qJD(3) * t93;
t126 = qJD(3) * t87;
t134 = t87 * t88;
t99 = (t88 * t126 + (-t93 * t94 + t134) * qJD(2)) * pkin(1);
t24 = qJD(1) * t99 - t125 * t62;
t132 = t88 * t93;
t100 = t87 * t94 + t132;
t98 = (qJD(2) * t100 + t125 * t88) * pkin(1);
t25 = qJD(1) * t98 + t126 * t62;
t5 = qJD(9) * t14 - t24 * t89 - t25 * t83;
t138 = t49 * t89;
t15 = t44 * t83 - t138;
t6 = qJD(9) * t15 + t24 * t83 - t25 * t89;
t159 = -t5 * t83 - t6 * t89;
t61 = pkin(3) * t80 + t74;
t86 = sin(qJ(4));
t92 = cos(qJ(4));
t48 = t116 * t92 + t61 * t86;
t81 = sin(qJ(10));
t141 = t48 * t81;
t46 = -t116 * t86 + t61 * t92;
t76 = qJD(4) + t80;
t43 = pkin(4) * t76 + t46;
t82 = cos(qJ(10));
t12 = -t43 * t82 + t141;
t123 = qJD(4) * t92;
t135 = t86 * t88;
t101 = t92 * t94 - t135;
t124 = qJD(4) * t86;
t97 = (qJD(2) * t101 - t124 * t88) * pkin(1);
t22 = qJD(1) * t97 + t123 * t61;
t133 = t88 * t92;
t102 = -t86 * t94 - t133;
t96 = (qJD(2) * t102 - t123 * t88) * pkin(1);
t23 = qJD(1) * t96 - t124 * t61;
t3 = qJD(10) * t12 - t22 * t82 - t23 * t81;
t140 = t48 * t82;
t13 = t43 * t81 + t140;
t4 = qJD(10) * t13 + t22 * t81 - t23 * t82;
t158 = -t3 * t81 - t4 * t82;
t157 = pkin(1) * (qJD(2) + qJD(6));
t156 = m(4) / 0.2e1;
t155 = m(5) / 0.2e1;
t1 = t4 * mrSges(11,1);
t150 = mrSges(5,1) * t23 + t1;
t2 = t6 * mrSges(10,1);
t149 = mrSges(4,1) * t25 + t2;
t147 = mrSges(4,1) * t77;
t146 = mrSges(5,1) * t76;
t69 = qJD(9) + t77;
t137 = t69 * mrSges(10,1);
t131 = t162 * qJD(1);
t68 = qJD(10) + t76;
t127 = t68 * mrSges(11,1);
t122 = qJD(9) * t83;
t121 = qJD(9) * t89;
t70 = pkin(3) * t92 + pkin(4);
t119 = qJD(10) * t70;
t118 = qJD(10) * t86;
t41 = t103 * t157;
t31 = qJD(1) * t41;
t42 = t104 * t157;
t32 = qJD(1) * t42;
t115 = mrSges(7,1) * t32 - t31 * mrSges(7,2);
t78 = t94 * pkin(1);
t73 = t78 + pkin(2);
t114 = pkin(1) * t134 - t73 * t93;
t72 = t78 + pkin(3);
t109 = -pkin(1) * t135 + t72 * t92;
t50 = pkin(4) + t109;
t58 = pkin(1) * t133 + t72 * t86;
t108 = -t50 * t82 + t58 * t81;
t107 = t50 * t81 + t58 * t82;
t51 = pkin(6) + t114;
t59 = -pkin(1) * t132 - t73 * t87;
t106 = -t51 * t89 + t59 * t83;
t105 = t51 * t83 + t59 * t89;
t95 = -t3 * mrSges(11,2) - t24 * mrSges(4,2) - t22 * mrSges(5,2) - t5 * mrSges(10,2) + t115 + t149 + t150;
t79 = qJD(1) + qJD(8);
t71 = -pkin(2) * t93 + pkin(6);
t57 = -t74 * t93 + t64;
t56 = t101 * t130;
t54 = t100 * t130;
t53 = t102 * t130;
t40 = t126 * t73 + t98;
t39 = -t125 * t73 + t99;
t38 = -t124 * t72 + t96;
t37 = t123 * t72 + t97;
t36 = t71 * t122 + (-t87 * t121 + (-t83 * t93 - t87 * t89) * qJD(3)) * pkin(2);
t35 = -t71 * t121 + (-t87 * t122 + (-t83 * t87 + t89 * t93) * qJD(3)) * pkin(2);
t34 = t81 * t119 + (t82 * t118 + (t81 * t92 + t82 * t86) * qJD(4)) * pkin(3);
t33 = -t82 * t119 + (t81 * t118 + (t81 * t86 - t82 * t92) * qJD(4)) * pkin(3);
t29 = -t54 * t83 - t57 * t89;
t28 = -t54 * t89 + t57 * t83;
t27 = -t53 * t81 - t56 * t82;
t26 = -t53 * t82 + t56 * t81;
t19 = -t47 * t89 - t139;
t18 = t47 * t83 - t138;
t17 = -t46 * t82 + t141;
t16 = t46 * t81 + t140;
t10 = qJD(9) * t105 + t39 * t83 - t40 * t89;
t9 = qJD(9) * t106 - t39 * t89 - t40 * t83;
t8 = qJD(10) * t107 + t37 * t81 - t38 * t82;
t7 = qJD(10) * t108 - t37 * t82 - t38 * t81;
t11 = [t131 - t39 * t77 * mrSges(4,2) + t40 * t147 - t7 * t68 * mrSges(11,2) - t9 * t69 * mrSges(10,2) + t10 * t137 - t41 * t143 + t42 * t136 - t37 * t76 * mrSges(5,2) + t38 * t146 + m(4) * (t114 * t25 + t24 * t59 - t49 * t39 + t47 * t40) + m(5) * (t109 * t23 + t22 * t58 + t48 * t37 + t46 * t38) + m(10) * (t14 * t10 - t105 * t5 + t106 * t6 - t15 * t9) + m(11) * (-t107 * t3 + t108 * t4 + t12 * t8 - t13 * t7) + t95 + m(7) * (-t55 * t41 + t52 * t42 + (t103 * t32 - t104 * t31) * pkin(1)) + t8 * t127 + t162 * t79 - (t80 + qJD(1)) * t163 * pkin(1) * qJD(2); m(11) * (t12 * t34 - t13 * t33 + t158 * t70) + m(10) * (t14 * t36 - t15 * t35 + t159 * t71) + ((t27 - t33) * mrSges(11,2) + (-t26 + t34) * mrSges(11,1)) * t68 + ((t29 - t35) * mrSges(10,2) + (-t28 + t36) * mrSges(10,1)) * t69 + (-t53 * mrSges(5,1) + t56 * mrSges(5,2) + (-mrSges(5,1) * t86 - mrSges(5,2) * t92) * qJD(4) * pkin(3)) * t76 + (-t54 * mrSges(4,1) + t57 * mrSges(4,2) + (mrSges(4,1) * t87 + mrSges(4,2) * t93) * qJD(3) * pkin(2)) * t77 - m(5) * (t46 * t53 + t48 * t56) - m(4) * (t47 * t54 - t49 * t57) - m(11) * (t12 * t26 - t13 * t27) - m(10) * (t14 * t28 - t15 * t29) + t95 + 0.2e1 * ((qJD(3) * t49 - t25) * t93 * t156 + (m(10) * (t5 * t89 - t6 * t83) / 0.2e1 + (qJD(3) * t47 - t24) * t156) * t87) * pkin(2) + 0.2e1 * ((qJD(4) * t48 + t23) * t92 * t155 + ((-qJD(4) * t46 + t22) * t155 + m(11) * (-t3 * t82 + t4 * t81) / 0.2e1) * t86) * pkin(3) + t163 * t130 * (-qJD(2) + t80) + t161; -m(10) * (t14 * t18 - t15 * t19) - t18 * t137 - t49 * t147 + (t19 * t69 - t5) * mrSges(10,2) + (t47 * t77 - t24) * mrSges(4,2) + (m(10) * ((t14 * t83 + t15 * t89) * qJD(9) + t159) + (mrSges(10,1) * t83 + mrSges(10,2) * t89) * qJD(9) * t69) * pkin(6) + t149; -m(11) * (t12 * t16 - t13 * t17) + t48 * t146 - t16 * t127 + (t46 * t76 - t22) * mrSges(5,2) + (t17 * t68 - t3) * mrSges(11,2) + (m(11) * ((t12 * t81 + t13 * t82) * qJD(10) + t158) + (mrSges(11,1) * t81 + mrSges(11,2) * t82) * qJD(10) * t68) * pkin(4) + t150; 0; t115 + t161; 0; -t130 * t164 * t79 + t131; -t15 * t137 + t2 + (t14 * t69 - t5) * mrSges(10,2); -t13 * t127 + t1 + (t12 * t68 - t3) * mrSges(11,2); 0; 0;];
tauc = t11(:);
