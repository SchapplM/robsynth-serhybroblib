% Calculate minimal parameter regressor of coriolis joint torque vector for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% tauc_reg [6x38]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh2m2OL_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:06:37
% EndTime: 2020-06-30 18:07:04
% DurationCPUTime: 4.25s
% Computational Cost: add. (7759->318), mult. (20455->457), div. (0->0), fcn. (17486->10), ass. (0->193)
t141 = sin(qJ(6));
t146 = cos(qJ(6));
t142 = sin(qJ(5));
t144 = sin(qJ(3));
t148 = cos(qJ(3));
t149 = cos(qJ(2));
t217 = qJD(1) * t149;
t145 = sin(qJ(2));
t218 = qJD(1) * t145;
t113 = t144 * t218 - t148 * t217;
t115 = -t144 * t217 - t148 * t218;
t143 = sin(qJ(4));
t147 = cos(qJ(4));
t184 = -t113 * t147 + t143 * t115;
t251 = cos(qJ(5));
t259 = t143 * t113 + t115 * t147;
t52 = t142 * t259 + t184 * t251;
t294 = qJD(6) + t52;
t138 = qJD(2) + qJD(3);
t137 = qJD(4) + t138;
t129 = qJD(5) + t137;
t269 = t142 * t184 - t251 * t259;
t236 = t141 * t269;
t42 = t129 * t146 + t236;
t277 = t294 * t42;
t213 = qJD(6) * t141;
t233 = t146 * t269;
t44 = -t141 * t129 + t233;
t41 = t44 * t213;
t117 = t144 * t145 - t148 * t149;
t118 = t144 * t149 + t145 * t148;
t182 = t117 * t143 - t118 * t147;
t274 = t138 * qJD(1);
t153 = t182 * t274;
t174 = t259 * qJD(4);
t152 = t174 + t153;
t203 = qJD(5) * t251;
t214 = qJD(5) * t142;
t172 = t117 * qJD(3);
t95 = -t117 * qJD(2) - t172;
t156 = t95 * qJD(1);
t157 = t118 * t274;
t215 = qJD(4) * t147;
t216 = qJD(4) * t143;
t38 = -t113 * t215 + t115 * t216 - t143 * t157 + t147 * t156;
t13 = t142 * t152 + t184 * t203 + t214 * t259 + t251 * t38;
t6 = -t42 * qJD(6) + t146 * t13;
t212 = qJD(6) * t146;
t7 = -t129 * t213 + t13 * t141 + t212 * t269;
t299 = (-t6 + t277) * t146 + t41 + (t44 * t52 + t7) * t141;
t234 = t146 * t294;
t14 = qJD(5) * t269 + t142 * t38 - t251 * t152;
t238 = t141 * t14;
t297 = -t234 * t294 + t269 * t44 + t238;
t235 = t146 * t44;
t245 = t6 * t141;
t296 = -t235 * t294 - t245;
t295 = t269 * t52;
t240 = pkin(4) * qJD(2);
t208 = t148 * t240;
t119 = pkin(2) * t138 + t208;
t191 = qJD(3) * t208;
t209 = t144 * t240;
t192 = t143 * t209;
t221 = (qJD(3) + qJD(4)) * t192;
t69 = (qJD(4) * t119 + t191) * t147 - t221;
t225 = t144 * t147;
t181 = -t143 * t148 - t225;
t163 = (t181 * qJD(3) - t144 * t215) * pkin(4);
t70 = qJD(2) * t163 - t119 * t216;
t198 = -t142 * t69 + t251 * t70;
t100 = t119 * t143 + t147 * t209;
t206 = t251 * t100;
t99 = t147 * t119 - t192;
t94 = pkin(5) * t137 + t99;
t60 = t142 * t94 + t206;
t22 = -t60 * qJD(5) + t198;
t133 = -pkin(4) * t149 - pkin(1);
t121 = t133 * qJD(1);
t98 = t113 * pkin(2) + t121;
t57 = -pkin(5) * t184 + t98;
t262 = t57 * t269;
t160 = t22 - t262;
t293 = t269 ^ 2 - t52 ^ 2;
t292 = -t129 * t52 + t13;
t165 = t100 * t214 - t142 * t70 - t94 * t203 - t251 * t69;
t291 = -t57 * t52 + t165;
t45 = t294 * t213;
t200 = t14 * t146 + t45;
t290 = t141 * t294 * t52 - t269 * t42 + t200;
t27 = -pkin(3) * t52 + t57;
t18 = -t141 * t27 + t146 * t60;
t286 = t18 * t269;
t285 = t294 * t269;
t284 = pkin(3) * t269;
t283 = t129 * t269 - t14;
t186 = t141 * t60 + t146 * t27;
t282 = t22 * t146 - t186 * t269;
t252 = t259 * pkin(5);
t243 = t98 * t184;
t275 = -t69 - t243;
t26 = -t184 ^ 2 + t259 ^ 2;
t155 = t259 * t98 + t70;
t24 = -t137 * t259 + t152;
t210 = qJD(1) * qJD(2);
t265 = -0.2e1 * t210;
t244 = t259 * t184;
t23 = -t137 * t184 + t38;
t250 = pkin(4) * t145;
t92 = t147 * t117 + t118 * t143;
t55 = -t142 * t92 - t182 * t251;
t249 = t14 * t55;
t109 = t181 * t240;
t226 = t143 * t144;
t180 = t147 * t148 - t226;
t110 = t180 * t240;
t131 = pkin(2) * t147 + pkin(5);
t205 = t251 * t143;
t242 = -t251 * t109 + t142 * t110 - t131 * t214 + (-t143 * t203 + (-t142 * t147 - t205) * qJD(4)) * pkin(2);
t227 = t142 * t143;
t241 = t142 * t109 + t251 * t110 - t131 * t203 - (-t143 * t214 + (t251 * t147 - t227) * qJD(4)) * pkin(2);
t232 = t115 * t113;
t230 = t121 * t113;
t229 = t121 * t115;
t228 = t142 * t100;
t220 = t145 ^ 2 - t149 ^ 2;
t151 = qJD(1) ^ 2;
t219 = pkin(1) * t151;
t211 = t115 * pkin(2);
t136 = t145 * t240;
t207 = t294 * t212;
t202 = t145 * t210;
t128 = pkin(4) * t202;
t20 = -pkin(5) * t174 + t128 + (t118 * pkin(2) - pkin(5) * t182) * t274;
t5 = t14 * pkin(3) + t20;
t201 = -qJD(6) * t60 - t5;
t59 = t251 * t94 - t228;
t56 = pkin(3) * t129 + t59;
t197 = t294 * t56;
t135 = pkin(4) * t218;
t61 = -t211 - t252;
t29 = t61 + t284;
t132 = pkin(4) * t148 + pkin(2);
t107 = -pkin(4) * t226 + t132 * t147 + pkin(5);
t111 = pkin(4) * t225 + t132 * t143;
t74 = t142 * t107 + t251 * t111;
t196 = qJD(6) * t74 - t135 - t29;
t173 = t118 * qJD(3);
t96 = t118 * qJD(2) + t173;
t82 = t96 * pkin(2) + t136;
t108 = pkin(2) * t205 + t142 * t131;
t195 = qJD(6) * t108 - t29;
t193 = qJD(3) * (-qJD(2) - t138);
t190 = pkin(1) * t265;
t189 = t294 * t203;
t188 = (-qJD(3) + t138) * t240;
t185 = t143 * t95 + t147 * t96;
t164 = -t182 * qJD(4) + t185;
t178 = t142 * t182 - t251 * t92;
t40 = -t92 * qJD(4) - t143 * t96 + t147 * t95;
t15 = t178 * qJD(5) - t142 * t164 + t251 * t40;
t187 = -t15 * t294 + t249;
t102 = pkin(2) * t117 + t133;
t177 = (-t56 + t59) * t294;
t63 = t251 * t99 - t228;
t176 = (-t56 + t63) * t294;
t175 = t251 * t107 - t142 * t111;
t66 = pkin(5) * t92 + t102;
t31 = -pkin(3) * t178 + t66;
t170 = qJD(6) * t294 * t31 + t56 * t15 + t22 * t55;
t88 = t132 * t215 + (t180 * qJD(3) - t144 * t216) * pkin(4);
t89 = -t132 * t216 + t163;
t33 = t175 * qJD(5) + t142 * t89 + t251 * t88;
t169 = t14 * t74 - t294 * t33 - t197;
t168 = t108 * t14 + t241 * t294 - t197;
t16 = t55 * qJD(5) + t142 * t40 + t251 * t164;
t30 = pkin(5) * t164 + t82;
t162 = t56 * qJD(6) * t55 + t31 * t14 + t201 * t178 - (t16 * pkin(3) + t30) * t294;
t150 = qJD(2) ^ 2;
t130 = t251 * pkin(5) + pkin(3);
t104 = -pkin(2) * t227 + t251 * t131 + pkin(3);
t101 = t135 - t211;
t72 = pkin(3) + t175;
t71 = pkin(2) * t157 + t128;
t68 = -t113 ^ 2 + t115 ^ 2;
t65 = -t115 * t138 - t157;
t64 = t113 * t138 + t156;
t62 = -t142 * t99 - t206;
t58 = t135 + t61;
t34 = -t74 * qJD(5) - t142 * t88 + t251 * t89;
t32 = -t252 + t284;
t25 = t27 * t213;
t1 = [0, 0, 0, 0.2e1 * t149 * t202, t220 * t265, t150 * t149, -t150 * t145, 0, t145 * t190, t149 * t190, -t115 * t95 + t118 * t156, -t95 * t113 + t115 * t96 + (t117 ^ 2 - t118 ^ 2) * t274, t95 * t138, -t96 * t138, 0, t113 * t136 + t121 * t96 + (t133 * t173 + (t117 * t250 + t133 * t118) * qJD(2)) * qJD(1), -t115 * t136 + t121 * t95 + (-t133 * t172 + (-t133 * t117 + t118 * t250) * qJD(2)) * qJD(1), -t182 * t38 - t259 * t40, t40 * t184 + t259 * t185 - t38 * t92 + (-t153 - 0.2e1 * t174) * t182, t40 * t137, -t164 * t137, 0, -t82 * t184 + t71 * t92 + t98 * t185 + (-t102 * t259 - t182 * t98) * qJD(4) - t102 * t153, t102 * t38 - t182 * t71 - t259 * t82 + t40 * t98, t13 * t55 + t15 * t269, t13 * t178 + t15 * t52 - t16 * t269 - t249, t15 * t129, -t16 * t129, 0, t14 * t66 + t16 * t57 - t178 * t20 - t30 * t52, t13 * t66 + t15 * t57 + t20 * t55 + t269 * t30, -t55 * t41 + (t15 * t44 + t55 * t6) * t146, (-t141 * t44 - t146 * t42) * t15 + (-t245 - t146 * t7 + (t141 * t42 - t235) * qJD(6)) * t55, -t146 * t187 - t44 * t16 + t178 * t6 - t45 * t55, t141 * t187 + t42 * t16 - t178 * t7 - t207 * t55, -t14 * t178 - t16 * t294, t186 * t16 + t25 * t178 + (t165 * t178 + t170) * t141 + t162 * t146, t18 * t16 + (-(-qJD(6) * t27 - t165) * t178 + t170) * t146 - t162 * t141; 0, 0, 0, -t145 * t151 * t149, t220 * t151, 0, 0, 0, t145 * t219, t149 * t219, -t232, t68, t64, t65, 0, t229 + (-t113 * t218 + t144 * t193) * pkin(4), t230 + (t115 * t218 + t148 * t193) * pkin(4), t244, t26, t23, t24, 0, t101 * t184 + t89 * t137 + t155, t101 * t259 - t88 * t137 + t275, -t295, t293, t292, t283, 0, t34 * t129 + t52 * t58 + t160, -t33 * t129 - t269 * t58 + t291, t296, t299, t297, t290, t285, t141 * t169 - t196 * t234 + t34 * t42 + t72 * t7 + t282, -t286 + t34 * t44 + t72 * t6 + (t196 * t294 - t22) * t141 + t169 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232, t68, t64, t65, 0, t144 * t188 + t229, t148 * t188 + t230, t244, t26, t23, t24, 0, -t109 * t137 + (-t115 * t184 - t137 * t216) * pkin(2) + t155, -t259 * t211 + t110 * t137 - t243 + (-t191 + (-pkin(2) * t137 - t119) * qJD(4)) * t147 + t221, -t295, t293, t292, t283, 0, t242 * t129 + t52 * t61 + t160, t241 * t129 - t269 * t61 + t291, t296, t299, t297, t290, t285, t104 * t7 + t168 * t141 - t195 * t234 + t242 * t42 + t282, t104 * t6 - t286 + t242 * t44 + (t195 * t294 - t22) * t141 + t168 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, t26, t23, t24, 0, t100 * t137 + t155, t99 * t137 + t275, -t295, t293, t292, t283, 0, -t52 * t252 - t62 * t129 - t262 + (-t206 + (-pkin(5) * t129 - t94) * t142) * qJD(5) + t198, t63 * t129 + (-t129 * t203 + t259 * t269) * pkin(5) + t291, t296, t299, t297, t290, t285, t32 * t234 + t130 * t7 - t62 * t42 + t176 * t141 + (-t141 * t189 + (-qJD(5) * t42 - t207 + t238) * t142) * pkin(5) + t282, t130 * t6 - t286 - t62 * t44 + (-t294 * t32 - t22) * t141 + t176 * t146 + (-t146 * t189 + (-qJD(5) * t44 + t200) * t142) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t295, t293, t292, t283, 0, t60 * t129 + t160, t59 * t129 + t291, t296, t299, t297, t290, t285, t60 * t42 + (t233 * t294 + t7) * pkin(3) + t177 * t141 + t282, -t22 * t141 - t286 + t60 * t44 + (-t236 * t294 + t6) * pkin(3) + t177 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t42, -t42 ^ 2 + t44 ^ 2, t6 + t277, t294 * t44 - t7, -t14, t141 * t165 + t146 * t201 + t18 * t294 - t56 * t44 + t25, t141 * t5 + t146 * t165 + t56 * t42 + (qJD(6) - t294) * t186;];
tauc_reg = t1;
