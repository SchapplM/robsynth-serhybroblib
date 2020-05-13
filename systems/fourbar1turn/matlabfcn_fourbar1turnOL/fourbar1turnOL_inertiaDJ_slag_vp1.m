% Calculate time derivative of joint inertia matrix for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnOL_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:49
% EndTime: 2020-04-12 19:40:56
% DurationCPUTime: 3.95s
% Computational Cost: add. (2852->323), mult. (4766->518), div. (0->0), fcn. (3714->8), ass. (0->196)
t131 = qJD(2) + qJD(3);
t134 = qJ(2) + qJ(3);
t127 = cos(t134);
t126 = sin(t134);
t230 = Icges(4,4) * t126;
t162 = Icges(4,2) * t127 + t230;
t151 = t131 * t162;
t229 = Icges(4,4) * t127;
t168 = Icges(4,1) * t126 + t229;
t152 = t131 * t168;
t184 = rSges(4,1) * t126 + rSges(4,2) * t127;
t278 = t131 * t184;
t137 = sin(qJ(1));
t140 = cos(qJ(1));
t136 = sin(qJ(2));
t139 = cos(qJ(2));
t255 = rSges(3,1) * t139;
t186 = -rSges(3,2) * t136 + t255;
t279 = t140 * rSges(3,3) - t137 * t186;
t250 = rSges(4,2) * t126;
t253 = rSges(4,1) * t127;
t185 = t250 - t253;
t135 = sin(qJ(4));
t138 = cos(qJ(4));
t247 = rSges(5,2) * t138;
t115 = rSges(5,1) * t135 + t247;
t153 = qJD(4) * t115;
t277 = t137 * t153;
t231 = Icges(3,4) * t139;
t165 = -Icges(3,2) * t136 + t231;
t85 = Icges(3,6) * t137 + t140 * t165;
t232 = Icges(3,4) * t136;
t171 = Icges(3,1) * t139 - t232;
t89 = Icges(3,5) * t137 + t140 * t171;
t172 = t136 * t85 - t139 * t89;
t276 = t137 * t172;
t227 = Icges(5,4) * t138;
t161 = -Icges(5,2) * t135 + t227;
t83 = Icges(5,6) * t137 + t140 * t161;
t228 = Icges(5,4) * t135;
t167 = Icges(5,1) * t138 - t228;
t87 = Icges(5,5) * t137 + t140 * t167;
t176 = t135 * t83 - t138 * t87;
t275 = t137 * t176;
t163 = -Icges(4,2) * t126 + t229;
t67 = Icges(4,6) * t137 - t140 * t163;
t169 = Icges(4,1) * t127 - t230;
t69 = Icges(4,5) * t137 - t140 * t169;
t181 = t126 * t67 - t127 * t69;
t274 = t137 * t181;
t84 = -Icges(3,6) * t140 + t137 * t165;
t88 = -Icges(3,5) * t140 + t137 * t171;
t173 = t136 * t84 - t139 * t88;
t272 = t140 * t173;
t82 = -Icges(5,6) * t140 + t137 * t161;
t86 = -Icges(5,5) * t140 + t137 * t167;
t177 = t135 * t82 - t138 * t86;
t271 = t140 * t177;
t264 = Icges(4,5) * t140 + t137 * t169;
t265 = Icges(4,6) * t140 + t137 * t163;
t182 = -t126 * t265 + t127 * t264;
t270 = t140 * t182;
t204 = qJD(2) * t139;
t206 = qJD(1) * t140;
t269 = -t136 * t206 - t137 * t204;
t240 = t137 * rSges(3,3);
t93 = t140 * t186 + t240;
t268 = -t137 * t93 - t140 * t279;
t158 = Icges(4,5) * t127 - Icges(4,6) * t126;
t180 = -t126 * t162 + t127 * t168;
t267 = qJD(1) * t180 + t131 * t158;
t156 = Icges(5,5) * t138 - Icges(5,6) * t135;
t78 = -Icges(5,3) * t140 + t137 * t156;
t266 = Icges(4,3) * t140 + t137 * t158;
t159 = Icges(3,5) * t139 - Icges(3,6) * t136;
t80 = -Icges(3,3) * t140 + t137 * t159;
t263 = 2 * m(3);
t262 = 2 * m(4);
t261 = 2 * m(5);
t132 = t137 ^ 2;
t133 = t140 ^ 2;
t260 = t137 / 0.2e1;
t259 = -t140 / 0.2e1;
t116 = rSges(3,1) * t136 + rSges(3,2) * t139;
t258 = m(3) * t116;
t257 = m(5) * t115;
t256 = pkin(2) * t139;
t246 = rSges(4,3) * t140;
t72 = t137 * t185 - t246;
t129 = t137 * rSges(4,3);
t211 = t140 * t250 + t129;
t73 = -t140 * t253 + t211;
t34 = t137 * t72 + t140 * t73;
t252 = rSges(5,1) * t138;
t248 = rSges(5,2) * t135;
t245 = rSges(5,3) * t140;
t244 = t135 * t86;
t243 = t135 * t87;
t242 = t136 * t88;
t241 = t136 * t89;
t128 = t137 * rSges(5,3);
t238 = t138 * t82;
t237 = t138 * t83;
t236 = t139 * t84;
t235 = t139 * t85;
t65 = Icges(4,3) * t137 - t140 * t158;
t217 = qJD(1) * t65;
t79 = Icges(5,3) * t137 + t140 * t156;
t216 = qJD(1) * t79;
t81 = Icges(3,3) * t137 + t140 * t159;
t215 = qJD(1) * t81;
t214 = t126 * t131;
t213 = t127 * t131;
t207 = qJD(1) * t137;
t210 = rSges(5,3) * t206 + t207 * t248;
t209 = t140 * t252 + t128;
t208 = t132 + t133;
t205 = qJD(2) * t136;
t203 = qJD(2) * t140;
t202 = qJD(4) * t135;
t201 = qJD(4) * t138;
t200 = qJD(4) * t140;
t199 = t140 * t248;
t198 = pkin(2) * t205;
t195 = -pkin(2) * t136 + t184;
t37 = qJD(1) * t265 + t140 * t151;
t192 = t131 * t69 + t37;
t38 = qJD(1) * t67 + t137 * t151;
t191 = t131 * t264 - t38;
t39 = qJD(1) * t264 + t140 * t152;
t190 = t131 * t67 - t39;
t40 = qJD(1) * t69 + t137 * t152;
t189 = -t131 * t265 - t40;
t188 = t132 / 0.2e1 + t133 / 0.2e1;
t11 = t137 * t182 + t140 * t266;
t12 = -t140 * t65 + t274;
t13 = -t137 * t266 + t270;
t14 = t137 * t65 + t140 * t181;
t157 = Icges(4,5) * t126 + Icges(4,6) * t127;
t150 = t131 * t157;
t35 = qJD(1) * t266 + t140 * t150;
t36 = t137 * t150 + t217;
t187 = -t140 * ((t140 * t36 + (t12 - t270) * qJD(1)) * t140 + (t11 * qJD(1) + (t126 * t37 - t127 * t39 + t67 * t213 + t69 * t214 + t217) * t137 + (qJD(1) * t181 + t191 * t126 - t189 * t127 - t35) * t140) * t137) + t137 * ((t137 * t35 + (t13 - t274) * qJD(1)) * t137 + (t14 * qJD(1) + (-t126 * t38 + t127 * t40 + t213 * t265 + t214 * t264) * t140 + (-t36 + t190 * t127 + t192 * t126 + (t182 + t65) * qJD(1)) * t137) * t140) + (-t11 * t140 + t12 * t137) * t207 + (-t13 * t140 + t14 * t137) * t206;
t183 = -t248 + t252;
t170 = Icges(3,1) * t136 + t231;
t166 = Icges(5,1) * t135 + t227;
t164 = Icges(3,2) * t139 + t232;
t160 = Icges(5,2) * t138 + t228;
t155 = -pkin(1) - t183;
t142 = rSges(4,3) * t206 + t140 * t278 - t185 * t207;
t10 = -t73 * t207 + t137 * (t137 * t278 + (t140 * t185 + t129) * qJD(1)) + t140 * t142 + t72 * t206;
t75 = t163 * t131;
t76 = t169 * t131;
t141 = -qJD(1) * t157 + (t76 - t151) * t127 + (-t75 - t152) * t126;
t154 = (t126 * t190 - t127 * t192 - t267 * t137 + t140 * t141) * t260 + (t126 * t189 + t127 * t191 + t141 * t137 + t140 * t267) * t259 + (t126 * t264 + t127 * t265 + t137 * t180 + t140 * t157) * t207 / 0.2e1 + (-t126 * t69 - t127 * t67 - t137 * t157 + t140 * t180) * t206 / 0.2e1;
t149 = qJD(2) * t170;
t148 = qJD(2) * t164;
t147 = qJD(2) * (-Icges(3,5) * t136 - Icges(3,6) * t139);
t146 = qJD(4) * t166;
t145 = qJD(4) * t160;
t144 = qJD(4) * (-Icges(5,5) * t135 - Icges(5,6) * t138);
t143 = -t185 - t256;
t107 = t186 * qJD(2);
t106 = t183 * qJD(4);
t92 = -t199 + t209;
t90 = t137 * t183 - t245;
t77 = t185 * t131;
t71 = (pkin(1) - t248) * t140 + t209;
t70 = t137 * t155 + t245;
t63 = t195 * t140;
t62 = t195 * t137;
t58 = (-t253 + t256) * t140 + t211;
t57 = t137 * t143 + t246;
t56 = -t137 * rSges(3,1) * t205 + (t140 * t255 + t240) * qJD(1) + t269 * rSges(3,2);
t55 = qJD(1) * t279 - t116 * t203;
t46 = t137 * t147 + t215;
t45 = -qJD(1) * t80 + t140 * t147;
t44 = t137 * t144 + t216;
t43 = -qJD(1) * t78 + t140 * t144;
t42 = t277 + (t140 * t155 - t128) * qJD(1);
t41 = -t200 * t247 - pkin(1) * t207 + (-t135 * t200 - t138 * t207) * rSges(5,1) + t210;
t31 = pkin(2) * t269 - t137 * t77 + t184 * t206;
t30 = -t184 * t207 - t140 * t77 + (t136 * t207 - t139 * t203) * pkin(2);
t25 = t208 * t256 + t34;
t24 = (t198 - t278) * t137 + (t140 * t143 - t129) * qJD(1);
t23 = (-t136 * t203 - t139 * t207) * pkin(2) + t142;
t22 = t137 * t81 - t140 * t172;
t21 = t137 * t80 - t272;
t20 = t137 * t79 - t140 * t176;
t19 = t137 * t78 - t271;
t18 = -t140 * t81 - t276;
t17 = -t137 * t173 - t140 * t80;
t16 = -t140 * t79 - t275;
t15 = -t137 * t177 - t140 * t78;
t9 = -t198 * t208 + t10;
t1 = [(-t279 * t56 + t55 * t93) * t263 + t168 * t213 + t126 * t76 - t162 * t214 + t127 * t75 + (t23 * t58 + t24 * t57) * t262 + (t41 * t71 + t42 * t70) * t261 + (t171 - t164) * t205 + (t170 + t165) * t204 + (t167 - t160) * t202 + (t166 + t161) * t201; (-t172 * qJD(2) + t136 * (-qJD(1) * t88 - t140 * t149) + t139 * (-qJD(1) * t84 - t140 * t148)) * t260 + (-t173 * qJD(2) + t136 * (qJD(1) * t89 - t137 * t149) + t139 * (qJD(1) * t85 - t137 * t148)) * t259 + m(3) * ((-t137 * t55 + t140 * t56) * t116 + t268 * t107) + m(4) * (t23 * t62 + t24 * t63 + t30 * t57 + t31 * t58) + t188 * t159 * qJD(2) + ((t241 / 0.2e1 + t235 / 0.2e1 - t93 * t258) * t140 + (t242 / 0.2e1 + t236 / 0.2e1 + t279 * t258) * t137) * qJD(1) + t154; ((-t137 * t279 + t140 * t93) * (qJD(1) * t268 + t137 * t56 + t140 * t55) + t208 * t116 * t107) * t263 + (t22 * t137 - t140 * t21) * t206 + t137 * ((t137 * t45 + (t21 + t276) * qJD(1)) * t137 + (t22 * qJD(1) + (t204 * t84 + t205 * t88) * t140 + (-t46 + (-t235 - t241) * qJD(2) + (-t173 + t81) * qJD(1)) * t137) * t140) + (t18 * t137 - t140 * t17) * t207 - t140 * ((t140 * t46 + (t18 + t272) * qJD(1)) * t140 + (t17 * qJD(1) + (-t204 * t85 - t205 * t89 + t215) * t137 + (-t45 + (t236 + t242) * qJD(2) - t172 * qJD(1)) * t140) * t137) + (t25 * t9 + t30 * t63 + t31 * t62) * t262 + t187; m(4) * ((-t137 * t58 - t140 * t57) * t77 - (-t137 * t23 - t140 * t24 + (t137 * t57 - t140 * t58) * qJD(1)) * t184) + t154; m(4) * (t10 * t25 + t34 * t9 + (-t137 * t62 - t140 * t63) * t77 - (-t137 * t31 - t140 * t30 + (t137 * t63 - t140 * t62) * qJD(1)) * t184) + t187; (-t184 * t208 * t77 + t10 * t34) * t262 + t187; (-qJD(4) * t176 + t135 * (-qJD(1) * t86 - t140 * t146) + t138 * (-qJD(1) * t82 - t140 * t145)) * t260 + (-qJD(4) * t177 + t135 * (qJD(1) * t87 - t137 * t146) + t138 * (qJD(1) * t83 - t137 * t145)) * t259 + m(5) * ((-t137 * t41 - t140 * t42) * t115 + (-t137 * t71 - t140 * t70) * t106) + t188 * t156 * qJD(4) + ((t243 / 0.2e1 + t237 / 0.2e1 - t71 * t257) * t140 + (t244 / 0.2e1 + t238 / 0.2e1 + t70 * t257) * t137) * qJD(1); 0; 0; ((t137 * t90 + t140 * t92) * ((qJD(1) * t90 - t140 * t153 + t210) * t140 + (-t277 + (-t199 - t92 + t128) * qJD(1)) * t137) + t208 * t115 * t106) * t261 + (t20 * t137 - t140 * t19) * t206 + t137 * ((t137 * t43 + (t19 + t275) * qJD(1)) * t137 + (t20 * qJD(1) + (t201 * t82 + t202 * t86) * t140 + (-t44 + (-t237 - t243) * qJD(4) + (-t177 + t79) * qJD(1)) * t137) * t140) + (t16 * t137 - t140 * t15) * t207 - t140 * ((t140 * t44 + (t16 + t271) * qJD(1)) * t140 + (t15 * qJD(1) + (-t201 * t83 - t202 * t87 + t216) * t137 + (-t43 + (t238 + t244) * qJD(4) - t176 * qJD(1)) * t140) * t137); 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
