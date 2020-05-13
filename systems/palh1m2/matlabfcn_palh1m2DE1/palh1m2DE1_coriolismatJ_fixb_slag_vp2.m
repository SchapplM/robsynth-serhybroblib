% Calculate matrix of centrifugal and coriolis load on the joints for
% palh1m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh1m2DE1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE1_coriolismatJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:57:39
% EndTime: 2020-05-01 20:57:47
% DurationCPUTime: 1.81s
% Computational Cost: add. (2633->268), mult. (2714->330), div. (0->0), fcn. (1370->90), ass. (0->164)
t228 = (qJ(3) + qJ(2));
t204 = 2 * t228;
t272 = sin(t204) / 0.2e1;
t164 = 2 * qJ(2);
t271 = sin(t164) / 0.2e1;
t138 = pkin(22) + pkin(21);
t208 = (pkin(18) - t138);
t200 = (-pkin(20) + t208);
t109 = (qJ(4) + t200);
t97 = qJ(2) + t109;
t110 = (-qJ(4) + t200);
t98 = -qJ(2) + t110;
t270 = cos(t98) / 0.4e1 + cos(t97) / 0.4e1;
t142 = sin(pkin(20));
t144 = cos(pkin(20));
t152 = sin(pkin(18));
t156 = cos(pkin(18));
t154 = cos(qJ(3));
t155 = cos(qJ(2));
t229 = t154 * t155;
t115 = pkin(5) * t229;
t150 = sin(qJ(3));
t123 = pkin(5) * t150 + pkin(1);
t151 = sin(qJ(2));
t36 = t123 * t151 - t115;
t230 = t151 * t154;
t37 = pkin(5) * t230 + t123 * t155;
t15 = t152 * t37 - t156 * t36;
t211 = t152 * t36 + t37 * t156;
t269 = t142 * t15 + t144 * t211;
t70 = t150 * t155 + t230;
t240 = t152 * t70;
t231 = t150 * t151;
t69 = t229 - t231;
t210 = t69 * t156 + t240;
t239 = t156 * t70;
t22 = t152 * t69 - t239;
t268 = t142 * t22 + t144 * t210;
t267 = -t142 * t211 + t144 * t15;
t266 = -t142 * t210 + t144 * t22;
t193 = qJ(2) + t200;
t187 = qJ(3) + t193;
t182 = qJ(4) + t187;
t194 = -qJ(2) + t200;
t189 = -qJ(3) + t194;
t185 = -qJ(4) + t189;
t261 = -sin(t185) / 0.4e1 + sin(t182) / 0.4e1;
t260 = sin(t98) / 0.4e1 - sin(t97) / 0.4e1;
t136 = qJ(2) + pkin(17) - pkin(18);
t259 = 0.2e1 * t136;
t258 = -pkin(9) / 0.2e1;
t257 = m(5) + m(6);
t251 = pkin(2) * mrSges(10,1);
t250 = pkin(5) * mrSges(6,1);
t249 = pkin(5) * mrSges(6,2);
t248 = pkin(9) * mrSges(6,2);
t247 = cos(t136);
t139 = qJ(3) + pkin(19);
t246 = cos(t139);
t244 = pkin(9) * t152;
t243 = pkin(9) * t156;
t149 = sin(qJ(4));
t153 = cos(qJ(4));
t96 = mrSges(6,1) * t153 - mrSges(6,2) * t149;
t242 = t37 * t96;
t212 = t164 + t139;
t102 = pkin(2) * mrSges(10,2) * sin(t212);
t222 = t257 * pkin(5);
t114 = t222 + mrSges(4,1) + mrSges(11,1);
t120 = sin(t136);
t183 = -qJ(4) + t187;
t184 = qJ(4) + t189;
t169 = -sin(t183) / 0.4e1 + sin(t184) / 0.4e1;
t170 = cos(t184) / 0.4e1 + cos(t183) / 0.4e1;
t180 = cos(t182);
t181 = cos(t185);
t202 = -qJ(2) + t208;
t196 = -qJ(3) + t202;
t201 = qJ(2) + t208;
t197 = qJ(3) + t201;
t213 = qJ(2) + t139;
t198 = 0.2e1 * t213;
t206 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t218 = pkin(9) * m(6) + mrSges(5,1);
t166 = ((-cos(t196) / 0.2e1 - cos(t197) / 0.2e1) * mrSges(11,2) + (sin(t196) / 0.2e1 - sin(t197) / 0.2e1) * mrSges(11,1)) * pkin(4) + (t222 * t272 + (-t180 / 0.4e1 - t181 / 0.4e1 + t170) * mrSges(6,2) + (t169 - t261) * mrSges(6,1) + (sin(t189) / 0.2e1 - sin(t187) / 0.2e1) * t218 + (cos(t187) / 0.2e1 - cos(t189) / 0.2e1) * t206) * pkin(5) - Ifges(9,4) * cos(t198) - (Ifges(11,4) + Ifges(4,4)) * cos(t204) + (pkin(2) ^ 2 * m(10) - Ifges(9,1) + Ifges(9,2)) * sin(t198) / 0.2e1 + (-Ifges(4,1) + Ifges(4,2) - Ifges(11,1) + Ifges(11,2)) * t272;
t188 = qJ(2) + t110;
t190 = -qJ(2) + t109;
t171 = sin(t190) / 0.4e1 - sin(t188) / 0.4e1;
t172 = -cos(t188) / 0.4e1 - cos(t190) / 0.4e1;
t116 = sin(t213);
t117 = cos(t213);
t128 = sin(t228);
t131 = cos(t228);
t147 = mrSges(11,2) + mrSges(4,2);
t177 = mrSges(9,2) * t117 + t114 * t128 + (pkin(2) * m(10) + mrSges(9,1)) * t116 + t147 * t131;
t203 = cos(t212);
t223 = pkin(18) - pkin(22);
t214 = qJ(2) + t223;
t215 = -qJ(2) + t223;
t226 = t164 + qJ(3);
t216 = cos(t226);
t99 = pkin(1) * t147 * sin(t226);
t1 = ((mrSges(3,1) + mrSges(10,1)) * t155 + t177) * pkin(15) - t102 - t99 - pkin(15) * (mrSges(3,2) + mrSges(10,2)) * t151 + t203 * t251 + t166 + Ifges(7,4) * cos(t259) + (-Ifges(7,2) + Ifges(7,1)) * sin(t259) / 0.2e1 + (t114 * t216 + (sin(t214) / 0.2e1 + sin(t215) / 0.2e1) * mrSges(8,2) + (-cos(t214) / 0.2e1 - cos(t215) / 0.2e1) * mrSges(8,1) + (-cos(t201) / 0.2e1 - cos(t202) / 0.2e1) * m(11) * pkin(4) + (t171 - t260) * mrSges(6,2) + (t172 - t270) * mrSges(6,1) + (-cos(t193) / 0.2e1 - cos(t194) / 0.2e1) * t218 + (pkin(1) * t271 - pkin(15) * t155) * (-m(4) - m(8) - m(11) - t257) + (-sin(t193) / 0.2e1 - sin(t194) / 0.2e1) * t206) * pkin(1) - (-Ifges(3,4) - Ifges(10,4)) * cos(t164) + (Ifges(3,1) + Ifges(10,1) - Ifges(3,2) - Ifges(10,2)) * t271 + (-mrSges(7,1) * t247 + mrSges(7,2) * t120) * pkin(14);
t241 = t1 * qJD(1);
t205 = 2 * t200;
t111 = qJ(4) + t205;
t112 = -qJ(4) + t205;
t132 = pkin(11) * mrSges(6,2) - Ifges(6,6);
t133 = pkin(11) * mrSges(6,1) - Ifges(6,5);
t159 = pkin(9) * mrSges(6,1);
t163 = 0.2e1 * qJ(4);
t76 = 2 * t109;
t77 = 2 * t110;
t2 = (sin(t163) / 0.4e1 + sin(t77) / 0.8e1 - sin(t76) / 0.8e1) * (Ifges(6,2) - Ifges(6,1)) + (-cos(t163) / 0.2e1 + cos(t76) / 0.4e1 + cos(t77) / 0.4e1) * Ifges(6,4) + (t153 * t258 + (cos(t109) / 0.2e1 + cos(t110) / 0.2e1) * pkin(15) + t170 * pkin(5)) * mrSges(6,2) + (t133 - t248) * cos(t111) / 0.4e1 - (t133 + t248) * cos(t112) / 0.4e1 - (t159 + t132) * sin(t111) / 0.4e1 + (t159 - t132) * sin(t112) / 0.4e1 + (t149 * t258 + (sin(t109) / 0.2e1 - sin(t110) / 0.2e1) * pkin(15) + t169 * pkin(5)) * mrSges(6,1) + t261 * t250 + (t181 + t180) * t249 / 0.4e1 + ((t172 + t270) * mrSges(6,1) + (t171 + t260) * mrSges(6,2)) * pkin(1);
t238 = t2 * qJD(1);
t217 = mrSges(10,2) * sin(t139);
t232 = t147 * t150;
t3 = -t102 / 0.2e1 - t99 / 0.2e1 - pkin(2) * t217 / 0.2e1 + (t203 / 0.2e1 - t246 / 0.2e1) * t251 + (t232 / 0.2e1 + (-t154 / 0.2e1 + t216 / 0.2e1) * t114) * pkin(1) + t177 * pkin(15) + t166;
t237 = t3 * qJD(1);
t236 = qJD(4) * t96;
t140 = qJ(4) + qJ(2);
t126 = sin(t140);
t227 = qJ(4) - qJ(2);
t127 = -sin(t227);
t129 = cos(t140);
t130 = cos(t227);
t134 = qJ(3) + t140;
t118 = sin(t134);
t135 = qJ(3) - t227;
t119 = sin(t135);
t121 = cos(t134);
t122 = cos(t135);
t176 = (t122 / 0.4e1 - t121 / 0.4e1) * mrSges(6,2) + (-t119 / 0.4e1 - t118 / 0.4e1) * mrSges(6,1);
t173 = t176 * pkin(5);
t168 = t173 + ((-t127 / 0.4e1 + t126 / 0.4e1) * mrSges(6,2) + (-t130 / 0.4e1 - t129 / 0.4e1) * mrSges(6,1)) * pkin(1);
t10 = -t242 / 0.2e1 + t168;
t235 = t10 * qJD(1);
t11 = -Ifges(9,6) * t116 - t128 * (Ifges(11,6) + Ifges(4,6)) + (-pkin(2) * mrSges(10,3) + Ifges(9,5)) * t117 + (-pkin(5) * mrSges(5,3) + Ifges(11,5) + Ifges(4,5)) * t131 + (-t118 / 0.2e1 + t119 / 0.2e1) * t250 - (t121 + t122) * t249 / 0.2e1;
t234 = t11 * qJD(3);
t191 = t70 * t96;
t13 = (-t191 / 0.2e1 + t176) * pkin(5);
t233 = t13 * qJD(1);
t124 = sin(t138);
t125 = cos(t138);
t95 = mrSges(6,1) * t149 + mrSges(6,2) * t153;
t224 = (t124 * t266 + t268 * t125) * t95 * pkin(5);
t209 = -t224 / 0.2e1;
t195 = mrSges(6,1) * t244 - t133 * t156;
t192 = (t114 * t154 - t232) * pkin(1);
t18 = (mrSges(10,1) * t246 + t217) * pkin(2) + t192;
t174 = -t156 * t231 + t240;
t175 = t152 * t231 + t239;
t167 = t95 * (((t142 * t152 + t144 * t156) * t125 - (t142 * t156 - t144 * t152) * t124) * t115 + ((-t142 * t175 + t144 * t174) * t125 - (t142 * t174 + t144 * t175) * t124) * pkin(5));
t4 = t209 + t167 / 0.2e1;
t186 = t18 * qJD(2) - t4 * qJD(4);
t145 = cos(pkin(19));
t143 = sin(pkin(19));
t64 = -mrSges(6,2) * t244 + t132 * t156;
t63 = mrSges(6,1) * t243 + t133 * t152;
t62 = mrSges(6,2) * t243 + t132 * t152;
t12 = pkin(5) * t191 / 0.2e1 + t173;
t9 = t242 / 0.2e1 + t168;
t5 = -t167 / 0.2e1 + t209;
t6 = [-qJD(2) * t1 - qJD(3) * t3 + qJD(4) * t2, t9 * qJD(4) + t234 - t241 + (-Ifges(7,5) * t120 - Ifges(7,6) * t247 - (Ifges(3,5) + Ifges(10,5)) * t151 - t155 * (Ifges(3,6) + Ifges(10,6)) + t11 + (-(-mrSges(11,3) - mrSges(4,3) - mrSges(5,3) - mrSges(8,3)) * t151 + (t127 / 0.2e1 + t126 / 0.2e1) * mrSges(6,2) + (t130 / 0.2e1 - t129 / 0.2e1) * mrSges(6,1)) * pkin(1)) * qJD(2), t11 * qJD(2) + t12 * qJD(4) + t234 - t237, t238 + t9 * qJD(2) + t12 * qJD(3) + ((-(t142 * t195 + t63 * t144) * t149 + (t142 * t64 - t144 * t62) * t153) * t125 + (-(-t142 * t63 + t144 * t195) * t149 + (t142 * t62 + t144 * t64) * t153) * t124 - (-pkin(15) + t36) * t95) * qJD(4); qJD(4) * t10 + t241, t18 * qJD(3), (t192 + ((mrSges(10,1) * t145 + mrSges(10,2) * t143) * t154 - (mrSges(10,1) * t143 - mrSges(10,2) * t145) * t150) * pkin(2)) * qJD(3) + t186, t235 - t4 * qJD(3) - (t267 * t124 + t269 * t125) * t236; qJD(4) * t13 + t237, -t186, 0, t233 + t4 * qJD(2) - (t124 * t268 - t266 * t125) * pkin(5) * t236; -qJD(2) * t10 - qJD(3) * t13 - t238, -t235 - t95 * (-t269 * t124 + t267 * t125) * qJD(2) + t5 * qJD(3), qJD(2) * t5 - qJD(3) * t224 - t233, 0;];
Cq = t6;
