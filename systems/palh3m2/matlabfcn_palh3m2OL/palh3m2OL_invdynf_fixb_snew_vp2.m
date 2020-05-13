% Calculate vector of cutting forces with Newton-Euler
% palh3m2OL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% qJDD [10x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x9]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = palh3m2OL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_invdynf_fixb_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_invdynf_fixb_snew_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m2OL_invdynf_fixb_snew_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_invdynf_fixb_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_invdynf_fixb_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_invdynf_fixb_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2OL_invdynf_fixb_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:55
% EndTime: 2020-05-07 04:34:13
% DurationCPUTime: 2.87s
% Computational Cost: add. (28866->256), mult. (63491->349), div. (0->0), fcn. (48752->18), ass. (0->137)
t169 = sin(qJ(8));
t128 = sin(qJ(4));
t136 = cos(qJ(4));
t120 = qJDD(2) + qJDD(3);
t129 = sin(qJ(3));
t137 = cos(qJ(3));
t130 = sin(qJ(2));
t138 = cos(qJ(2));
t140 = qJD(1) ^ 2;
t131 = sin(qJ(1));
t139 = cos(qJ(1));
t154 = -g(1) * t139 - g(2) * t131;
t98 = -pkin(12) * t140 + t154;
t153 = -t138 * g(3) - t130 * t98;
t82 = (t130 * t138 * t140 + qJDD(2)) * pkin(1) + t153;
t159 = -t130 * g(3) + t138 * t98;
t83 = (-t138 ^ 2 * t140 - qJD(2) ^ 2) * pkin(1) + t159;
t155 = t129 * t83 - t137 * t82;
t91 = (t129 * t130 - t137 * t138) * qJD(1);
t93 = (-t129 * t138 - t130 * t137) * qJD(1);
t38 = (t91 * t93 + t120) * pkin(4) + t155;
t122 = qJD(2) + qJD(3);
t151 = -t129 * t82 - t137 * t83;
t45 = (-t122 ^ 2 - t91 ^ 2) * pkin(4) + t151;
t168 = t128 * t38 + t136 * t45;
t125 = sin(qJ(7));
t133 = cos(qJ(7));
t167 = t125 * t82 + t133 * t83;
t126 = sin(qJ(6));
t166 = qJD(1) * t126;
t165 = qJD(1) * t130;
t134 = cos(qJ(6));
t164 = qJD(1) * t134;
t163 = qJD(1) * t138;
t162 = qJD(1) * qJD(2);
t161 = qJD(1) * qJD(6);
t121 = qJD(2) + qJD(7);
t119 = qJDD(2) + qJDD(7);
t156 = -t125 * t83 + t133 * t82;
t113 = qJDD(8) + t119;
t115 = qJD(8) + t121;
t90 = (-t125 * t130 + t133 * t138) * qJD(1);
t92 = (t125 * t138 + t130 * t133) * qJD(1);
t123 = sin(pkin(15));
t124 = cos(pkin(15));
t132 = cos(qJ(8));
t94 = t123 * t132 - t124 * t169;
t95 = -t123 * t169 - t124 * t132;
t55 = t90 * t95 - t92 * t94;
t103 = qJDD(1) * t130 + t138 * t162;
t158 = t130 * t162;
t105 = qJDD(1) * t138 - t158;
t61 = -qJD(7) * t92 - t103 * t125 + t105 * t133;
t63 = qJD(7) * t90 + t103 * t133 + t105 * t125;
t27 = qJD(8) * t55 + t61 * t94 + t63 * t95;
t118 = t121 ^ 2;
t67 = (-t123 * t92 - t124 * t90) * pkin(3);
t33 = -t92 * t67 + (t118 * t123 + t119 * t124) * pkin(3) + t156;
t34 = t90 * t67 + (-t118 * t124 + t119 * t123) * pkin(3) + t167;
t56 = t90 * t94 + t92 * t95;
t36 = -mrSges(9,1) * t55 + mrSges(9,2) * t56;
t53 = -mrSges(9,2) * t115 + mrSges(9,3) * t55;
t17 = m(9) * (t33 * t95 - t34 * t94) - t27 * mrSges(9,3) + t113 * mrSges(9,1) - t56 * t36 + t115 * t53;
t26 = -qJD(8) * t56 + t61 * t95 - t63 * t94;
t54 = mrSges(9,1) * t115 - mrSges(9,3) * t56;
t18 = m(9) * (t33 * t94 + t34 * t95) + t26 * mrSges(9,3) - t113 * mrSges(9,2) + t55 * t36 - t115 * t54;
t75 = -mrSges(8,1) * t90 + mrSges(8,2) * t92;
t84 = -mrSges(8,2) * t121 + mrSges(8,3) * t90;
t10 = m(8) * t156 + mrSges(8,1) * t119 - mrSges(8,3) * t63 + t121 * t84 + t17 * t95 + t18 * t94 - t75 * t92;
t101 = (-mrSges(3,1) * t138 + mrSges(3,2) * t130) * qJD(1);
t109 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t163;
t86 = mrSges(8,1) * t121 - mrSges(8,3) * t92;
t11 = m(8) * t167 - mrSges(8,2) * t119 + mrSges(8,3) * t61 - t121 * t86 - t17 * t94 + t18 * t95 + t75 * t90;
t114 = qJDD(4) + t120;
t116 = qJD(4) + t122;
t112 = t116 ^ 2;
t152 = -t128 * t45 + t136 * t38;
t127 = sin(qJ(5));
t135 = cos(qJ(5));
t62 = -qJD(3) * t93 + t103 * t129 - t105 * t137;
t64 = qJD(3) * t91 - t103 * t137 - t105 * t129;
t71 = -t128 * t93 + t136 * t91;
t32 = qJD(4) * t71 + t128 * t62 + t136 * t64;
t72 = t128 * t91 + t136 * t93;
t58 = t116 * t127 + t135 * t72;
t22 = -qJD(5) * t58 + t114 * t135 - t127 * t32;
t57 = t116 * t135 - t127 * t72;
t23 = qJD(5) * t57 + t114 * t127 + t135 * t32;
t68 = qJD(5) - t71;
t42 = -mrSges(6,2) * t68 + mrSges(6,3) * t57;
t43 = mrSges(6,1) * t68 - mrSges(6,3) * t58;
t48 = -pkin(8) * t71 - pkin(10) * t72;
t145 = m(6) * (-pkin(8) * t114 - pkin(10) * t112 + t48 * t72 - t152) - t22 * mrSges(6,1) + t23 * mrSges(6,2) - t57 * t42 + t58 * t43;
t47 = -mrSges(5,1) * t71 + mrSges(5,2) * t72;
t65 = -mrSges(5,2) * t116 + mrSges(5,3) * t71;
t12 = m(5) * t152 + mrSges(5,1) * t114 - mrSges(5,3) * t32 + t116 * t65 - t47 * t72 - t145;
t76 = -mrSges(4,1) * t91 + mrSges(4,2) * t93;
t85 = -mrSges(4,2) * t122 + mrSges(4,3) * t91;
t157 = t131 * g(1) - t139 * g(2);
t96 = -qJDD(1) * pkin(12) - t157;
t79 = t96 + (-t105 + t158) * pkin(1);
t144 = t79 + (t122 * t93 - t62) * pkin(4);
t31 = -qJD(4) * t72 - t128 * t64 + t136 * t62;
t19 = (-t116 * t71 - t32) * pkin(10) + (t116 * t72 - t31) * pkin(8) + t144;
t21 = -pkin(8) * t112 + pkin(10) * t114 + t48 * t71 + t168;
t30 = qJDD(5) - t31;
t39 = -mrSges(6,1) * t57 + mrSges(6,2) * t58;
t15 = m(6) * (-t127 * t21 + t135 * t19) - t23 * mrSges(6,3) + t30 * mrSges(6,1) - t58 * t39 + t68 * t42;
t16 = m(6) * (t127 * t19 + t135 * t21) + t22 * mrSges(6,3) - t30 * mrSges(6,2) + t57 * t39 - t68 * t43;
t66 = mrSges(5,1) * t116 - mrSges(5,3) * t72;
t9 = m(5) * t168 - mrSges(5,2) * t114 + mrSges(5,3) * t31 - t116 * t66 - t127 * t15 + t135 * t16 + t47 * t71;
t6 = m(4) * t155 + mrSges(4,1) * t120 - mrSges(4,3) * t64 + t12 * t136 + t122 * t85 + t128 * t9 - t76 * t93;
t87 = mrSges(4,1) * t122 - mrSges(4,3) * t93;
t7 = m(4) * t151 - mrSges(4,2) * t120 + mrSges(4,3) * t62 - t12 * t128 - t122 * t87 + t136 * t9 + t76 * t91;
t4 = m(3) * t153 + qJDD(2) * mrSges(3,1) - t103 * mrSges(3,3) + qJD(2) * t109 + t133 * t10 - t101 * t165 + t125 * t11 - t129 * t7 - t137 * t6;
t107 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t165;
t5 = m(3) * t159 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t105 - qJD(2) * t107 - t10 * t125 + t101 * t163 + t11 * t133 + t129 * t6 - t137 * t7;
t100 = (-mrSges(7,1) * t134 + mrSges(7,2) * t126) * qJD(1);
t102 = qJDD(1) * t126 + t134 * t161;
t108 = -qJD(6) * mrSges(7,2) + mrSges(7,3) * t164;
t99 = pkin(6) * t140 + t154;
t51 = m(7) * (-g(3) * t134 - t126 * t99) - t102 * mrSges(7,3) + qJDD(6) * mrSges(7,1) - t100 * t166 + qJD(6) * t108;
t104 = qJDD(1) * t134 - t126 * t161;
t106 = qJD(6) * mrSges(7,1) - mrSges(7,3) * t166;
t52 = m(7) * (-g(3) * t126 + t134 * t99) + t104 * mrSges(7,3) - qJDD(6) * mrSges(7,2) + t100 * t164 - qJD(6) * t106;
t160 = t126 * t52 + t130 * t5 + t134 * t51 + t138 * t4;
t150 = t106 * t126 - t108 * t134;
t149 = t107 * t130 - t109 * t138;
t148 = m(7) * (qJDD(1) * pkin(6) - t157) - t104 * mrSges(7,1) + t102 * mrSges(7,2);
t147 = t26 * mrSges(9,1) - t27 * mrSges(9,2) + t55 * t53 - m(9) * ((-t123 * t63 - t124 * t61 + (-t123 * t90 + t124 * t92) * t121) * pkin(3) + t79) - t56 * t54;
t146 = -m(5) * t144 + t31 * mrSges(5,1) - mrSges(5,2) * t32 - t127 * t16 - t135 * t15 + t71 * t65 - t72 * t66;
t143 = -m(8) * t79 + t61 * mrSges(8,1) - t63 * mrSges(8,2) + t90 * t84 - t92 * t86 + t147;
t142 = -m(4) * t79 + t62 * mrSges(4,1) - t64 * mrSges(4,2) + t91 * t85 - t93 * t87 + t146;
t141 = m(3) * t96 - t105 * mrSges(3,1) + t103 * mrSges(3,2) - t142 - t143;
t8 = (-t149 - t150) * qJD(1) + qJDD(1) * mrSges(2,1) - t140 * mrSges(2,2) + m(2) * t157 - t141 - t148;
t1 = m(2) * t154 - mrSges(2,1) * t140 - qJDD(1) * mrSges(2,2) - t126 * t51 - t130 * t4 + t134 * t52 + t138 * t5;
t2 = [-m(1) * g(1) + t1 * t139 - t131 * t8, t1, t5, t7, t9, t16, t52, t11, t18, 0, 0, 0, 0; -m(1) * g(2) + t1 * t131 + t139 * t8, t8, t4, t6, t12, t15, t51, t10, t17, 0, 0, 0, 0; (-m(1) - m(2)) * g(3) + t160, -m(2) * g(3) + t160, qJD(1) * t149 + t141, -t142, -t146, t145, qJD(1) * t150 + t148, -t143, -t147, 0, 0, 0, 0;];
f_new = t2;
