% Calculate vector of cutting forces with Newton-Euler
% palh1m1OL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% qJDD [13x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% f_new [3x11]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = palh1m1OL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_invdynf_fixb_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_invdynf_fixb_snew_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m1OL_invdynf_fixb_snew_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_invdynf_fixb_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_invdynf_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_invdynf_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_invdynf_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:48
% EndTime: 2020-04-15 19:29:07
% DurationCPUTime: 3.58s
% Computational Cost: add. (33501->324), mult. (73487->436), div. (0->0), fcn. (56398->22), ass. (0->169)
t168 = sin(qJ(4));
t177 = cos(qJ(4));
t169 = sin(qJ(3));
t170 = sin(qJ(2));
t178 = cos(qJ(3));
t179 = cos(qJ(2));
t121 = (t169 * t179 + t170 * t178) * qJD(1);
t124 = (-t169 * t170 + t178 * t179) * qJD(1);
t156 = qJDD(2) + qJDD(3);
t181 = qJD(1) ^ 2;
t171 = sin(qJ(1));
t180 = cos(qJ(1));
t196 = -t180 * g(1) - t171 * g(2);
t130 = -t181 * pkin(15) + t196;
t117 = t170 * g(3) - t179 * t130;
t105 = (-t170 * t179 * t181 + qJDD(2)) * pkin(1) + t117;
t116 = -t179 * g(3) - t170 * t130;
t106 = (-t170 ^ 2 * t181 - qJD(2) ^ 2) * pkin(1) + t116;
t209 = t169 * t105 + t178 * t106;
t44 = (t121 * t124 + t156) * pkin(5) + t209;
t159 = qJD(2) + qJD(3);
t198 = -t178 * t105 + t169 * t106;
t51 = (-t124 ^ 2 - t159 ^ 2) * pkin(5) + t198;
t212 = t168 * t44 + t177 * t51;
t211 = sin(qJ(10));
t165 = sin(qJ(7));
t174 = cos(qJ(7));
t210 = t165 * t105 + t174 * t106;
t164 = sin(qJ(8));
t173 = cos(qJ(8));
t208 = t173 * t116 + t164 * t117;
t166 = sin(qJ(6));
t207 = qJD(1) * t166;
t206 = qJD(1) * t170;
t175 = cos(qJ(6));
t205 = qJD(1) * t175;
t204 = qJD(1) * t179;
t203 = qJD(1) * qJD(2);
t202 = qJD(1) * qJD(6);
t158 = qJD(2) + qJD(7);
t157 = qJD(2) + qJD(8);
t155 = qJDD(2) + qJDD(7);
t154 = qJDD(2) + qJDD(8);
t201 = t179 * t203;
t200 = t171 * g(1) - t180 * g(2);
t199 = t174 * t105 - t165 * t106;
t197 = -t164 * t116 + t173 * t117;
t123 = (-t165 * t170 + t174 * t179) * qJD(1);
t114 = t158 * mrSges(8,1) - t123 * mrSges(8,3);
t120 = (-t165 * t179 - t170 * t174) * qJD(1);
t160 = sin(pkin(19));
t161 = cos(pkin(19));
t162 = cos(qJ(10));
t126 = t160 * t162 - t161 * t211;
t127 = -t160 * t211 - t161 * t162;
t144 = qJDD(10) + t155;
t148 = qJD(10) + t158;
t67 = t127 * t120 - t126 * t123;
t135 = -t170 * qJDD(1) - t201;
t137 = t179 * qJDD(1) - t170 * t203;
t75 = -t123 * qJD(7) + t174 * t135 - t165 * t137;
t78 = t120 * qJD(7) + t165 * t135 + t174 * t137;
t29 = t67 * qJD(10) + t126 * t75 + t127 * t78;
t153 = t158 ^ 2;
t84 = (-t120 * t161 - t123 * t160) * pkin(4);
t39 = -t123 * t84 + (t153 * t160 + t155 * t161) * pkin(4) + t199;
t40 = t120 * t84 + (-t153 * t161 + t155 * t160) * pkin(4) + t210;
t68 = t126 * t120 + t127 * t123;
t42 = -t67 * mrSges(11,1) + t68 * mrSges(11,2);
t63 = -t148 * mrSges(11,2) + t67 * mrSges(11,3);
t18 = m(11) * (-t126 * t40 + t127 * t39) - t29 * mrSges(11,3) + t144 * mrSges(11,1) - t68 * t42 + t148 * t63;
t28 = -t68 * qJD(10) - t126 * t78 + t127 * t75;
t64 = t148 * mrSges(11,1) - t68 * mrSges(11,3);
t19 = m(11) * (t126 * t39 + t127 * t40) + t28 * mrSges(11,3) - t144 * mrSges(11,2) + t67 * t42 - t148 * t64;
t98 = -t120 * mrSges(8,1) + t123 * mrSges(8,2);
t10 = m(8) * t210 - t155 * mrSges(8,2) + t75 * mrSges(8,3) - t158 * t114 + t120 * t98 - t126 * t18 + t127 * t19;
t133 = (mrSges(3,1) * t170 + mrSges(3,2) * t179) * qJD(1);
t141 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t204;
t119 = (-t164 * t179 - t170 * t173) * qJD(1);
t110 = -t157 * mrSges(9,2) + t119 * mrSges(9,3);
t122 = (-t164 * t170 + t173 * t179) * qJD(1);
t163 = sin(qJ(9));
t172 = cos(qJ(9));
t146 = qJDD(9) + t154;
t149 = qJD(9) + t157;
t74 = -t122 * qJD(8) + t173 * t135 - t164 * t137;
t77 = t119 * qJD(8) + t164 * t135 + t173 * t137;
t89 = -t172 * t119 + t163 * t122;
t37 = t89 * qJD(9) - t163 * t74 - t172 * t77;
t52 = (t119 * t122 + t154) * pkin(2) + t197;
t91 = -t163 * t119 - t172 * t122;
t57 = -t89 * mrSges(10,1) + t91 * mrSges(10,2);
t62 = (-t119 ^ 2 - t157 ^ 2) * pkin(2) + t208;
t80 = -t149 * mrSges(10,2) + t89 * mrSges(10,3);
t20 = m(10) * (t163 * t62 - t172 * t52) - t37 * mrSges(10,3) + t146 * mrSges(10,1) - t91 * t57 + t149 * t80;
t35 = -t91 * qJD(9) + t163 * t77 - t172 * t74;
t82 = t149 * mrSges(10,1) - t91 * mrSges(10,3);
t21 = m(10) * (-t163 * t52 - t172 * t62) + t35 * mrSges(10,3) - t146 * mrSges(10,2) + t89 * t57 - t149 * t82;
t97 = -t119 * mrSges(9,1) + t122 * mrSges(9,2);
t16 = m(9) * t197 + t154 * mrSges(9,1) - t77 * mrSges(9,3) + t157 * t110 - t122 * t97 - t163 * t21 - t172 * t20;
t113 = t157 * mrSges(9,1) - t122 * mrSges(9,3);
t17 = m(9) * t208 - t154 * mrSges(9,2) + t74 * mrSges(9,3) - t157 * t113 + t119 * t97 + t163 * t20 - t172 * t21;
t147 = qJDD(4) + t156;
t150 = qJD(4) + t159;
t145 = t150 ^ 2;
t194 = -t168 * t51 + t177 * t44;
t167 = sin(qJ(5));
t176 = cos(qJ(5));
t76 = t124 * qJD(3) - t178 * t135 + t169 * t137;
t79 = -t121 * qJD(3) + t169 * t135 + t178 * t137;
t90 = -t168 * t121 + t177 * t124;
t38 = t90 * qJD(4) + t168 * t79 + t177 * t76;
t92 = t177 * t121 + t168 * t124;
t70 = t167 * t150 + t176 * t92;
t25 = -t70 * qJD(5) + t176 * t147 - t167 * t38;
t69 = t176 * t150 - t167 * t92;
t26 = t69 * qJD(5) + t167 * t147 + t176 * t38;
t85 = qJD(5) - t90;
t48 = -t85 * mrSges(6,2) + t69 * mrSges(6,3);
t49 = t85 * mrSges(6,1) - t70 * mrSges(6,3);
t59 = -t90 * pkin(9) - t92 * pkin(11);
t187 = m(6) * (-t147 * pkin(9) - t145 * pkin(11) + t92 * t59 - t194) - t25 * mrSges(6,1) + t26 * mrSges(6,2) - t69 * t48 + t70 * t49;
t58 = -t90 * mrSges(5,1) + t92 * mrSges(5,2);
t81 = -t150 * mrSges(5,2) + t90 * mrSges(5,3);
t11 = m(5) * t194 + t147 * mrSges(5,1) - t38 * mrSges(5,3) + t150 * t81 - t92 * t58 - t187;
t112 = t159 * mrSges(4,1) - t121 * mrSges(4,3);
t128 = -qJDD(1) * pkin(15) - t200;
t101 = t128 + (-t135 + t201) * pkin(1);
t186 = t101 + (t121 * t159 - t79) * pkin(5);
t36 = -t92 * qJD(4) - t168 * t76 + t177 * t79;
t22 = (-t150 * t90 - t38) * pkin(11) + (t150 * t92 - t36) * pkin(9) + t186;
t24 = -t145 * pkin(9) + t147 * pkin(11) + t90 * t59 + t212;
t34 = qJDD(5) - t36;
t45 = -t69 * mrSges(6,1) + t70 * mrSges(6,2);
t14 = m(6) * (-t167 * t24 + t176 * t22) - t26 * mrSges(6,3) + t34 * mrSges(6,1) - t70 * t45 + t85 * t48;
t15 = m(6) * (t167 * t22 + t176 * t24) + t25 * mrSges(6,3) - t34 * mrSges(6,2) + t69 * t45 - t85 * t49;
t83 = t150 * mrSges(5,1) - t92 * mrSges(5,3);
t8 = m(5) * t212 - t147 * mrSges(5,2) + t36 * mrSges(5,3) - t167 * t14 + t176 * t15 - t150 * t83 + t90 * t58;
t96 = -t124 * mrSges(4,1) + t121 * mrSges(4,2);
t6 = m(4) * t198 - t156 * mrSges(4,2) + t79 * mrSges(4,3) - t168 * t11 - t159 * t112 + t124 * t96 + t177 * t8;
t115 = -t159 * mrSges(4,2) + t124 * mrSges(4,3);
t7 = m(4) * t209 + t156 * mrSges(4,1) - t76 * mrSges(4,3) + t177 * t11 + t159 * t115 - t121 * t96 + t168 * t8;
t111 = -t158 * mrSges(8,2) + t120 * mrSges(8,3);
t9 = m(8) * t199 + t155 * mrSges(8,1) - t78 * mrSges(8,3) + t158 * t111 - t123 * t98 + t126 * t19 + t127 * t18;
t3 = m(3) * t116 - qJDD(2) * mrSges(3,2) + t135 * mrSges(3,3) - qJD(2) * t141 + t174 * t10 - t133 * t206 - t164 * t16 - t165 * t9 + t169 * t6 + t173 * t17 + t178 * t7;
t139 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t206;
t4 = m(3) * t117 + qJDD(2) * mrSges(3,1) - t137 * mrSges(3,3) + qJD(2) * t139 + t165 * t10 - t133 * t204 + t173 * t16 + t164 * t17 + t169 * t7 + t174 * t9 - t178 * t6;
t131 = t181 * pkin(14) + t196;
t132 = (-mrSges(7,1) * t175 + mrSges(7,2) * t166) * qJD(1);
t134 = t166 * qJDD(1) + t175 * t202;
t140 = -qJD(6) * mrSges(7,2) + mrSges(7,3) * t205;
t65 = m(7) * (-t175 * g(3) - t166 * t131) - t134 * mrSges(7,3) + qJDD(6) * mrSges(7,1) - t132 * t207 + qJD(6) * t140;
t136 = t175 * qJDD(1) - t166 * t202;
t138 = qJD(6) * mrSges(7,1) - mrSges(7,3) * t207;
t66 = m(7) * (-t166 * g(3) + t175 * t131) + t136 * mrSges(7,3) - qJDD(6) * mrSges(7,2) + t132 * t205 - qJD(6) * t138;
t195 = t166 * t66 - t170 * t4 + t175 * t65 + t179 * t3;
t193 = t166 * t138 - t175 * t140;
t192 = t170 * t139 + t179 * t141;
t191 = m(7) * (qJDD(1) * pkin(14) - t200) - t136 * mrSges(7,1) + t134 * mrSges(7,2);
t190 = t35 * mrSges(10,1) - t37 * mrSges(10,2) + t89 * t80 - t91 * t82 - m(10) * ((t122 * t157 - t74) * pkin(2) + t128);
t189 = t28 * mrSges(11,1) - t29 * mrSges(11,2) + t67 * t63 - m(11) * ((-t160 * t78 - t161 * t75 + (-t120 * t160 + t123 * t161) * t158) * pkin(4) + t101) - t68 * t64;
t188 = -m(5) * t186 + t36 * mrSges(5,1) - t38 * mrSges(5,2) - t176 * t14 - t167 * t15 + t90 * t81 - t92 * t83;
t185 = -m(8) * t101 + t75 * mrSges(8,1) - t78 * mrSges(8,2) + t120 * t111 - t123 * t114 + t189;
t184 = -m(9) * t128 + t74 * mrSges(9,1) - t77 * mrSges(9,2) + t119 * t110 - t122 * t113 + t190;
t183 = -m(4) * t101 + t79 * mrSges(4,1) - t76 * mrSges(4,2) - t121 * t112 + t124 * t115 + t188;
t182 = m(3) * t128 - t135 * mrSges(3,1) + t137 * mrSges(3,2) - t183 - t184 - t185;
t5 = (-t192 - t193) * qJD(1) - t181 * mrSges(2,2) + m(2) * t200 + qJDD(1) * mrSges(2,1) - t182 - t191;
t1 = m(2) * t196 - t181 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t166 * t65 - t170 * t3 + t175 * t66 - t179 * t4;
t2 = [-m(1) * g(1) + t180 * t1 - t171 * t5, t1, t3, t6, t8, t15, t66, t10, t17, t21, t19, 0, 0, 0, 0, 0, 0; -m(1) * g(2) + t171 * t1 + t180 * t5, t5, t4, t7, t11, t14, t65, t9, t16, t20, t18, 0, 0, 0, 0, 0, 0; (-m(1) - m(2)) * g(3) + t195, -m(2) * g(3) + t195, t192 * qJD(1) + t182, -t183, -t188, t187, t193 * qJD(1) + t191, -t185, -t184, -t190, -t189, 0, 0, 0, 0, 0, 0;];
f_new = t2;
