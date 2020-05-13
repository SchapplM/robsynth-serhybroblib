% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [12x12]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = picker2Dm2OL_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp2: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:58
% EndTime: 2020-05-09 23:19:05
% DurationCPUTime: 1.55s
% Computational Cost: add. (3698->234), mult. (8270->327), div. (0->0), fcn. (6404->14), ass. (0->177)
t141 = sin(qJ(10));
t142 = cos(qJ(10));
t122 = (t141 * mrSges(11,1) + t142 * mrSges(11,2)) * pkin(4);
t143 = sin(qJ(9));
t148 = cos(qJ(9));
t124 = (t143 * mrSges(10,1) + t148 * mrSges(10,2)) * pkin(6);
t145 = sin(qJ(4));
t152 = cos(qJ(2));
t147 = sin(qJ(2));
t150 = cos(qJ(4));
t195 = t147 * t150;
t114 = (-t145 * t152 - t195) * pkin(1);
t107 = t114 * mrSges(5,1);
t146 = sin(qJ(3));
t151 = cos(qJ(3));
t194 = t147 * t151;
t115 = (t146 * t152 + t194) * pkin(1);
t108 = t115 * mrSges(4,1);
t144 = sin(qJ(6));
t149 = cos(qJ(6));
t175 = (-(t144 * t147 - t149 * t152) * mrSges(7,2) + (t144 * t152 + t147 * t149) * mrSges(7,1)) * pkin(1);
t197 = t145 * t147;
t117 = (t150 * t152 - t197) * pkin(1);
t73 = -t141 * t114 - t142 * t117;
t210 = t73 * mrSges(11,2);
t130 = t146 * t147 * pkin(1);
t140 = t152 * pkin(1);
t118 = -t151 * t140 + t130;
t222 = t118 * mrSges(4,2);
t223 = t117 * mrSges(5,2);
t74 = -t143 * t115 - t148 * t118;
t237 = t74 * mrSges(10,2);
t167 = -t142 * t114 + t141 * t117;
t71 = t167 * mrSges(11,1);
t166 = -t148 * t115 + t143 * t118;
t72 = t166 * mrSges(10,1);
t259 = -(t147 * mrSges(3,1) + t152 * mrSges(3,2)) * pkin(1) + t107 + t108 + t175 - t210 - t222 - t223 - t237 + t71 + t72;
t75 = t175 * qJD(6);
t200 = t142 * t145;
t109 = (t141 * t150 + t200) * pkin(3);
t105 = t109 * mrSges(11,1);
t129 = t141 * t145 * pkin(3);
t241 = t150 * pkin(3);
t110 = -t142 * t241 + t129;
t205 = t110 * mrSges(11,2);
t216 = t150 * mrSges(5,2);
t220 = t145 * mrSges(5,1);
t257 = (t216 + t220) * pkin(3) - t105 + t205;
t176 = t151 * pkin(2) - pkin(6);
t196 = t146 * t148;
t102 = pkin(2) * t196 + t143 * t176;
t111 = (-t143 * t151 - t196) * pkin(2);
t198 = t143 * t146;
t116 = (t148 * t151 - t198) * pkin(2);
t224 = t116 * mrSges(10,2);
t225 = t111 * mrSges(10,1);
t164 = t225 / 0.2e1 - t224 / 0.2e1;
t215 = t151 * mrSges(4,2);
t219 = t146 * mrSges(4,1);
t139 = t140 + pkin(2);
t101 = -t151 * t139 + t130;
t227 = t101 * mrSges(4,2);
t104 = pkin(1) * t194 + t146 * t139;
t199 = t143 * t104;
t66 = -t148 * t101 - t199;
t238 = t66 * mrSges(10,2);
t193 = t148 * t104;
t65 = t143 * t101 - t193;
t58 = t65 * mrSges(10,1);
t234 = t58 / 0.2e1 - t238 / 0.2e1;
t250 = m(10) / 0.2e1;
t90 = pkin(6) + t101;
t61 = -t148 * t90 - t199;
t62 = t143 * t90 - t193;
t86 = t104 * mrSges(4,1);
t99 = -pkin(2) * t198 + t148 * t176;
t256 = t86 / 0.2e1 + (t102 * t66 + t111 * t61 - t116 * t62 + t99 * t65) * t250 - t227 / 0.2e1 + t164 + t234 + (t215 / 0.2e1 + t219 / 0.2e1) * pkin(2);
t242 = t105 / 0.2e1;
t173 = t242 - t205 / 0.2e1;
t138 = t140 + pkin(3);
t100 = -pkin(1) * t197 + t150 * t138;
t228 = t100 * mrSges(5,2);
t103 = pkin(1) * t195 + t145 * t138;
t203 = t141 * t103;
t64 = -t142 * t100 + t203;
t211 = t64 * mrSges(11,2);
t201 = t142 * t103;
t63 = t141 * t100 + t201;
t57 = t63 * mrSges(11,1);
t235 = t57 / 0.2e1 - t211 / 0.2e1;
t248 = m(11) / 0.2e1;
t89 = pkin(4) + t100;
t59 = -t142 * t89 + t203;
t60 = t141 * t89 + t201;
t85 = t103 * mrSges(5,1);
t137 = pkin(4) + t241;
t97 = -t142 * t137 + t129;
t98 = pkin(3) * t200 + t141 * t137;
t255 = -t85 / 0.2e1 + (t109 * t59 - t110 * t60 + t97 * t63 - t98 * t64) * t248 - t228 / 0.2e1 + t173 + t235 + (-t216 / 0.2e1 - t220 / 0.2e1) * pkin(3);
t254 = -pkin(4) / 0.2e1;
t252 = -pkin(6) / 0.2e1;
t249 = pkin(6) * m(10);
t247 = -t59 / 0.2e1;
t246 = -t61 / 0.2e1;
t245 = -t71 / 0.2e1;
t244 = -t72 / 0.2e1;
t243 = pkin(4) * m(11);
t240 = t61 * mrSges(10,2);
t239 = t62 * mrSges(10,1);
t236 = t99 * mrSges(10,2);
t233 = t71 / 0.2e1 - t210 / 0.2e1;
t232 = t72 / 0.2e1 - t237 / 0.2e1;
t208 = t98 * mrSges(11,1);
t209 = t97 * mrSges(11,2);
t231 = -t209 / 0.2e1 + t208 / 0.2e1;
t226 = t102 * mrSges(10,1);
t230 = -t236 / 0.2e1 - t226 / 0.2e1;
t1 = -m(5) * (t100 * t114 + t103 * t117) - m(4) * (t101 * t115 - t104 * t118) - m(11) * (t59 * t167 - t60 * t73) - m(10) * (t61 * t166 - t62 * t74) - t259;
t229 = t1 * qJD(1);
t213 = t59 * mrSges(11,2);
t212 = t60 * mrSges(11,1);
t161 = t211 - t57 + t85 + t228;
t10 = -m(11) * (t59 * t63 - t60 * t64) + t161;
t207 = t10 * qJD(1);
t162 = -t227 + t58 + t86 - t238;
t11 = -m(10) * (t61 * t65 - t62 * t66) - t162;
t206 = t11 * qJD(1);
t16 = t212 - t213;
t192 = t16 * qJD(1);
t17 = t239 - t240;
t191 = t17 * qJD(1);
t190 = t17 * qJD(9);
t33 = -t226 - t236;
t189 = t33 * qJD(9);
t188 = t175 * qJD(1);
t187 = t122 / 0.2e1;
t186 = t124 / 0.2e1;
t123 = (mrSges(9,1) * sin(qJ(8)) + mrSges(9,2) * cos(qJ(8))) * pkin(1);
t185 = t123 * qJD(1);
t184 = t16 * qJD(10);
t31 = t208 - t209;
t183 = t31 * qJD(10);
t182 = t239 / 0.2e1;
t177 = t212 / 0.2e1;
t155 = t107 / 0.2e1 - t223 / 0.2e1 + (-t141 * t73 - t142 * t167) * t243 / 0.2e1;
t2 = t210 / 0.2e1 + t245 - t155 + t255;
t30 = -m(11) * (t97 * t109 - t98 * t110) + t257;
t171 = t2 * qJD(1) - t30 * qJD(2);
t6 = t177 + t245 + (t247 + t73 / 0.2e1) * mrSges(11,2) + t231;
t170 = t6 * qJD(1) + t31 * qJD(2);
t158 = t225 + (t215 + t219) * pkin(2) - t224;
t32 = m(10) * (t102 * t116 + t99 * t111) + t158;
t156 = t108 / 0.2e1 - t222 / 0.2e1 + (-t143 * t74 - t148 * t166) * t249 / 0.2e1;
t5 = t237 / 0.2e1 + t244 - t156 + t256;
t169 = t5 * qJD(1) + t32 * qJD(2);
t8 = t182 + t244 + (t246 + t74 / 0.2e1) * mrSges(10,2) + t230;
t168 = t8 * qJD(1) + t33 * qJD(2);
t165 = -t240 / 0.2e1 + t182;
t163 = -t213 / 0.2e1 + t177;
t14 = t182 - t58 / 0.2e1 + (t246 + t66 / 0.2e1) * mrSges(10,2) + t186;
t28 = (-t116 / 0.2e1 + t99 / 0.2e1 + t148 * t252) * mrSges(10,2) + (t111 / 0.2e1 + t102 / 0.2e1 + t143 * t252) * mrSges(10,1);
t160 = -t14 * qJD(1) + t28 * qJD(2) - t124 * qJD(3);
t12 = t177 - t57 / 0.2e1 + (t247 + t64 / 0.2e1) * mrSges(11,2) + t187;
t26 = t242 + (-t98 / 0.2e1 + t141 * t254) * mrSges(11,1) + (-t110 / 0.2e1 + t97 / 0.2e1 + t142 * t254) * mrSges(11,2);
t159 = -t12 * qJD(1) + t26 * qJD(2) - t122 * qJD(4);
t121 = t124 * qJD(9);
t120 = t123 * qJD(8);
t119 = t122 * qJD(10);
t29 = t164 + t186 + t230;
t27 = t173 + t187 + t231;
t15 = t165 + t186 + t234;
t13 = t163 + t187 + t235;
t9 = t165 + t230 + t232;
t7 = t163 + t231 + t233;
t4 = t156 + t232 + t256;
t3 = t155 + t233 + t255;
t18 = [-t1 * qJD(2) - t11 * qJD(3) - t10 * qJD(4) + t120 + t184 + t190 + t75, t7 * qJD(10) + t4 * qJD(3) + t3 * qJD(4) + t9 * qJD(9) - t229 + t75 + (0.2e1 * (t102 * t74 + t99 * t166) * t250 + 0.2e1 * (t97 * t167 - t98 * t73) * t248 + m(5) * (t114 * t150 + t117 * t145) * pkin(3) + m(4) * (-t115 * t151 - t118 * t146) * pkin(2) + t259) * qJD(2), -t206 + t4 * qJD(2) + ((-t143 * t66 - t148 * t65) * t249 + t162) * qJD(3) + t15 * qJD(9), -t207 + t3 * qJD(2) + ((-t141 * t64 - t142 * t63) * t243 - t161) * qJD(4) + t13 * qJD(10), 0, qJD(2) * t175 + t188 + t75, 0, t120 + t185, t9 * qJD(2) + t15 * qJD(3) + t190 + t191, t7 * qJD(2) + t13 * qJD(4) + t184 + t192, 0, 0; t6 * qJD(10) + t5 * qJD(3) + t2 * qJD(4) + t8 * qJD(9) + t229, t32 * qJD(3) - t30 * qJD(4) + t183 + t189, ((-t111 * t148 - t116 * t143) * t249 + t158) * qJD(3) + t29 * qJD(9) + t169, ((-t109 * t142 - t110 * t141) * t243 - t257) * qJD(4) + t27 * qJD(10) + t171, 0, 0, 0, 0, t29 * qJD(3) + t168 + t189, t27 * qJD(4) + t170 + t183, 0, 0; -t5 * qJD(2) + t14 * qJD(9) + t206, -t28 * qJD(9) - t169, t121, 0, 0, 0, 0, 0, t121 - t160, 0, 0, 0; t12 * qJD(10) - t2 * qJD(2) + t207, -t26 * qJD(10) - t171, 0, t119, 0, 0, 0, 0, 0, t119 - t159, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t188, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t185, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t8 * qJD(2) - t14 * qJD(3) - t191, t28 * qJD(3) - t168, t160, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t6 * qJD(2) - t12 * qJD(4) - t192, t26 * qJD(4) - t170, 0, t159, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
Cq = t18;
