% Calculate matrix of centrifugal and coriolis load on the joints for
% picker2Dm1OL
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
% Cq [12x12]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = picker2Dm1OL_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_coriolismatJ_fixb_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_coriolismatJ_fixb_slag_vp2: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_coriolismatJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1OL_coriolismatJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1OL_coriolismatJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:45:01
% EndTime: 2020-05-11 05:45:06
% DurationCPUTime: 1.45s
% Computational Cost: add. (3698->238), mult. (8270->331), div. (0->0), fcn. (6404->14), ass. (0->179)
t140 = sin(qJ(10));
t141 = cos(qJ(10));
t121 = (t140 * mrSges(11,1) + t141 * mrSges(11,2)) * pkin(4);
t142 = sin(qJ(9));
t147 = cos(qJ(9));
t123 = (t142 * mrSges(10,1) + t147 * mrSges(10,2)) * pkin(6);
t144 = sin(qJ(4));
t151 = cos(qJ(2));
t146 = sin(qJ(2));
t149 = cos(qJ(4));
t190 = t146 * t149;
t113 = (-t144 * t151 - t190) * pkin(1);
t106 = t113 * mrSges(5,1);
t145 = sin(qJ(3));
t150 = cos(qJ(3));
t189 = t146 * t150;
t114 = (t145 * t151 + t189) * pkin(1);
t107 = t114 * mrSges(4,1);
t143 = sin(qJ(6));
t148 = cos(qJ(6));
t172 = (-(t143 * t146 - t148 * t151) * mrSges(7,2) + (t143 * t151 + t146 * t148) * mrSges(7,1)) * pkin(1);
t192 = t144 * t146;
t116 = (t149 * t151 - t192) * pkin(1);
t73 = -t140 * t113 - t141 * t116;
t205 = t73 * mrSges(11,2);
t129 = t145 * t146 * pkin(1);
t139 = t151 * pkin(1);
t117 = -t150 * t139 + t129;
t217 = t117 * mrSges(4,2);
t218 = t116 * mrSges(5,2);
t74 = -t142 * t114 - t147 * t117;
t234 = t74 * mrSges(10,2);
t72 = -t141 * t113 + t140 * t116;
t70 = t72 * mrSges(11,1);
t163 = -t147 * t114 + t142 * t117;
t71 = t163 * mrSges(10,1);
t257 = -(t146 * mrSges(3,1) + t151 * mrSges(3,2)) * pkin(1) + t106 + t107 + t172 - t205 - t217 - t218 - t234 + t70 + t71;
t75 = t172 * qJD(6);
t195 = t141 * t144;
t108 = (t140 * t149 + t195) * pkin(3);
t104 = t108 * mrSges(11,1);
t128 = t140 * t144 * pkin(3);
t238 = t149 * pkin(3);
t109 = -t141 * t238 + t128;
t201 = t109 * mrSges(11,2);
t225 = t104 / 0.2e1 - t201 / 0.2e1;
t211 = t149 * mrSges(5,2);
t215 = t144 * mrSges(5,1);
t255 = (t211 + t215) * pkin(3) - t104 + t201;
t173 = t150 * pkin(2) - pkin(6);
t191 = t145 * t147;
t101 = pkin(2) * t191 + t142 * t173;
t110 = (-t142 * t150 - t191) * pkin(2);
t193 = t142 * t145;
t115 = (t147 * t150 - t193) * pkin(2);
t219 = t115 * mrSges(10,2);
t220 = t110 * mrSges(10,1);
t162 = t220 / 0.2e1 - t219 / 0.2e1;
t210 = t150 * mrSges(4,2);
t214 = t145 * mrSges(4,1);
t138 = t139 + pkin(2);
t103 = pkin(1) * t189 + t145 * t138;
t221 = t103 * mrSges(4,1);
t100 = -t150 * t138 + t129;
t223 = t100 * mrSges(4,2);
t194 = t142 * t103;
t65 = -t147 * t100 - t194;
t235 = t65 * mrSges(10,2);
t188 = t147 * t103;
t64 = t142 * t100 - t188;
t236 = t64 * mrSges(10,1);
t230 = t236 / 0.2e1 - t235 / 0.2e1;
t249 = m(10) / 0.2e1;
t89 = pkin(6) + t100;
t60 = -t147 * t89 - t194;
t61 = t142 * t89 - t188;
t98 = -pkin(2) * t193 + t147 * t173;
t254 = (t101 * t65 + t110 * t60 - t115 * t61 + t98 * t64) * t249 - t223 / 0.2e1 + t221 / 0.2e1 + t162 + t230 + (t214 / 0.2e1 + t210 / 0.2e1) * pkin(2);
t253 = -pkin(4) / 0.2e1;
t251 = -pkin(6) / 0.2e1;
t248 = pkin(6) * m(10);
t247 = m(11) / 0.2e1;
t56 = t61 * mrSges(10,1);
t246 = t56 / 0.2e1;
t137 = t139 + pkin(3);
t102 = pkin(1) * t190 + t144 * t137;
t196 = t141 * t102;
t99 = -pkin(1) * t192 + t149 * t137;
t62 = t140 * t99 + t196;
t57 = t62 * mrSges(11,1);
t245 = t57 / 0.2e1;
t198 = t140 * t102;
t88 = pkin(4) + t99;
t58 = -t141 * t88 + t198;
t244 = t58 / 0.2e1;
t59 = t140 * t88 + t196;
t243 = -t59 / 0.2e1;
t242 = -t60 / 0.2e1;
t63 = -t141 * t99 + t198;
t241 = -t63 / 0.2e1;
t240 = -t71 / 0.2e1;
t239 = pkin(4) * m(11);
t237 = t60 * mrSges(10,2);
t233 = t98 * mrSges(10,2);
t232 = t99 * mrSges(5,2);
t207 = t59 * mrSges(11,1);
t208 = t58 * mrSges(11,2);
t231 = -t208 / 0.2e1 + t207 / 0.2e1;
t229 = t70 / 0.2e1 - t205 / 0.2e1;
t228 = t71 / 0.2e1 - t234 / 0.2e1;
t136 = pkin(4) + t238;
t97 = pkin(3) * t195 + t140 * t136;
t203 = t97 * mrSges(11,1);
t96 = -t141 * t136 + t128;
t204 = t96 * mrSges(11,2);
t227 = -t204 / 0.2e1 + t203 / 0.2e1;
t222 = t101 * mrSges(10,1);
t226 = -t233 / 0.2e1 - t222 / 0.2e1;
t1 = -m(4) * (t100 * t114 - t103 * t117) - m(11) * (t58 * t72 - t59 * t73) - m(10) * (t60 * t163 - t61 * t74) - m(5) * (t102 * t116 + t99 * t113) - t257;
t224 = t1 * qJD(1);
t206 = t63 * mrSges(11,2);
t85 = t102 * mrSges(5,1);
t161 = -t206 + t57 - t85 - t232;
t10 = -m(11) * (t58 * t62 - t59 * t63) - t161;
t202 = t10 * qJD(1);
t158 = t221 - t223 - t235 + t236;
t11 = m(10) * (t60 * t64 - t61 * t65) + t158;
t200 = t11 * qJD(1);
t16 = t207 - t208;
t187 = t16 * qJD(1);
t17 = -t56 + t237;
t186 = t17 * qJD(1);
t185 = t17 * qJD(9);
t33 = -t222 - t233;
t184 = t33 * qJD(9);
t183 = t172 * qJD(1);
t182 = t121 / 0.2e1;
t181 = t123 / 0.2e1;
t122 = (mrSges(9,1) * sin(qJ(8)) + mrSges(9,2) * cos(qJ(8))) * pkin(1);
t180 = t122 * qJD(1);
t179 = t16 * qJD(10);
t31 = t203 - t204;
t178 = t31 * qJD(10);
t177 = -t215 / 0.2e1;
t174 = -t206 / 0.2e1;
t170 = t182 + t227;
t169 = t246 - t237 / 0.2e1;
t155 = t245 - t85 / 0.2e1 + (t108 * t58 - t109 * t59 + t96 * t62 - t97 * t63) * t247 + t225;
t157 = -t106 / 0.2e1 - (-t140 * t73 - t141 * t72) * t239 / 0.2e1;
t2 = pkin(3) * t177 - t70 / 0.2e1 + (t241 + t73 / 0.2e1) * mrSges(11,2) + (-t99 / 0.2e1 - t238 / 0.2e1 + t116 / 0.2e1) * mrSges(5,2) + t155 + t157;
t30 = -m(11) * (t96 * t108 - t97 * t109) + t255;
t167 = t2 * qJD(1) - t30 * qJD(2);
t6 = (t244 + t96 / 0.2e1) * mrSges(11,2) + (t243 - t97 / 0.2e1) * mrSges(11,1) + t229;
t166 = -t6 * qJD(1) + t31 * qJD(2);
t156 = t220 + (t210 + t214) * pkin(2) - t219;
t32 = m(10) * (t101 * t115 + t98 * t110) + t156;
t153 = t107 / 0.2e1 - t217 / 0.2e1 + (-t142 * t74 - t147 * t163) * t248 / 0.2e1;
t4 = t234 / 0.2e1 + t240 - t153 + t254;
t165 = t4 * qJD(1) + t32 * qJD(2);
t8 = t246 + t240 + (t242 + t74 / 0.2e1) * mrSges(10,2) + t226;
t164 = t8 * qJD(1) + t33 * qJD(2);
t14 = t246 - t236 / 0.2e1 + (t242 + t65 / 0.2e1) * mrSges(10,2) + t181;
t28 = (-t115 / 0.2e1 + t98 / 0.2e1 + t147 * t251) * mrSges(10,2) + (t110 / 0.2e1 + t101 / 0.2e1 + t142 * t251) * mrSges(10,1);
t160 = -t14 * qJD(1) + t28 * qJD(2) - t123 * qJD(3);
t12 = t245 + (t140 * t253 + t243) * mrSges(11,1) + (t141 * t253 + t241 + t244) * mrSges(11,2);
t26 = t170 - t225;
t159 = t12 * qJD(1) - t26 * qJD(2) - t121 * qJD(4);
t120 = t123 * qJD(9);
t119 = t122 * qJD(8);
t118 = t121 * qJD(10);
t29 = t162 + t181 + t226;
t27 = t170 + t225;
t15 = t169 + t181 + t230;
t13 = t245 + t174 + t182 + t231;
t9 = t169 + t226 + t228;
t7 = t227 + t229 + t231;
t5 = t153 + t228 + t254;
t3 = -t232 / 0.2e1 + t174 - t218 / 0.2e1 + (t177 - t211 / 0.2e1) * pkin(3) + t155 - t157 + t229;
t18 = [-t1 * qJD(2) + t11 * qJD(3) - t10 * qJD(4) + t119 + t179 - t185 + t75, t7 * qJD(10) + t5 * qJD(3) + t3 * qJD(4) + t9 * qJD(9) - t224 + t75 + (0.2e1 * (t96 * t72 - t97 * t73) * t247 + 0.2e1 * (t101 * t74 + t163 * t98) * t249 + m(5) * (t113 * t149 + t116 * t144) * pkin(3) + m(4) * (-t114 * t150 - t117 * t145) * pkin(2) + t257) * qJD(2), t200 + t5 * qJD(2) + ((-t142 * t65 - t147 * t64) * t248 + t158) * qJD(3) + t15 * qJD(9), -t202 + t3 * qJD(2) + ((-t140 * t63 - t141 * t62) * t239 + t161) * qJD(4) + t13 * qJD(10), 0, qJD(2) * t172 + t183 + t75, 0, t119 + t180, t9 * qJD(2) + t15 * qJD(3) - t185 - t186, t7 * qJD(2) + t13 * qJD(4) + t179 + t187, 0, 0; -t6 * qJD(10) + t4 * qJD(3) + t2 * qJD(4) + t8 * qJD(9) + t224, t32 * qJD(3) - t30 * qJD(4) + t178 + t184, ((-t110 * t147 - t115 * t142) * t248 + t156) * qJD(3) + t29 * qJD(9) + t165, ((-t108 * t141 - t109 * t140) * t239 - t255) * qJD(4) + t27 * qJD(10) + t167, 0, 0, 0, 0, t29 * qJD(3) + t164 + t184, t27 * qJD(4) + t166 + t178, 0, 0; -t4 * qJD(2) + t14 * qJD(9) - t200, -t28 * qJD(9) - t165, t120, 0, 0, 0, 0, 0, t120 - t160, 0, 0, 0; -t12 * qJD(10) - t2 * qJD(2) + t202, t26 * qJD(10) - t167, 0, t118, 0, 0, 0, 0, 0, t118 - t159, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t183, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t180, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t8 * qJD(2) - t14 * qJD(3) + t186, t28 * qJD(3) - t164, t160, 0, 0, 0, 0, 0, 0, 0, 0, 0; t6 * qJD(2) + t12 * qJD(4) - t187, -t26 * qJD(4) - t166, 0, t159, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
Cq = t18;
