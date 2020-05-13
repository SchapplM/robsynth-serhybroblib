% Jacobian of explicit kinematic constraints of
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% W [12x4]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = palh3m1TE_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:16:33
% EndTime: 2020-04-17 15:16:37
% DurationCPUTime: 4.47s
% Computational Cost: add. (74765->140), mult. (112987->286), div. (4878->27), fcn. (71640->18), ass. (0->169)
t174 = sin(qJ(2));
t175 = sin(pkin(16));
t261 = cos(qJ(2));
t262 = cos(pkin(16));
t163 = t174 * t175 - t261 * t262;
t257 = pkin(5) * t163;
t224 = pkin(1) * t257;
t161 = -0.2e1 * t224;
t185 = pkin(5) ^ 2;
t190 = pkin(1) ^ 2;
t229 = t185 + t190;
t158 = t161 + t229;
t156 = 0.1e1 / t158;
t189 = 0.1e1 / pkin(2);
t234 = t156 * t189;
t183 = pkin(6) ^ 2;
t188 = pkin(2) ^ 2;
t155 = -t183 + t188 + t158;
t159 = pkin(1) - t257;
t230 = t161 + t190;
t273 = -pkin(6) - pkin(2);
t150 = (pkin(5) - t273) * (pkin(5) + t273) + t230;
t272 = -pkin(6) + pkin(2);
t151 = (pkin(5) - t272) * (pkin(5) + t272) + t230;
t238 = t151 * t150;
t192 = sqrt(-t238);
t164 = t174 * t262 + t175 * t261;
t256 = pkin(5) * t164;
t144 = t155 * t256 + t159 * t192;
t173 = sin(qJ(3));
t241 = t144 * t173;
t232 = t164 * t192;
t223 = pkin(5) * t232;
t143 = t155 * t159 - t223;
t177 = cos(qJ(3));
t242 = t143 * t177;
t137 = (-t242 / 0.2e1 + t241 / 0.2e1) * t234;
t240 = t144 * t177;
t243 = t143 * t173;
t138 = (t240 / 0.2e1 + t243 / 0.2e1) * t234;
t168 = pkin(18) + pkin(19);
t166 = sin(t168);
t167 = cos(t168);
t123 = t137 * t167 + t138 * t166;
t258 = t123 * pkin(4);
t225 = pkin(3) * t258;
t121 = 0.2e1 * t225;
t186 = pkin(4) ^ 2;
t231 = t121 + t186;
t271 = -pkin(8) - pkin(10);
t109 = (pkin(3) - t271) * (pkin(3) + t271) + t231;
t270 = -pkin(8) + pkin(10);
t110 = (pkin(3) - t270) * (pkin(3) + t270) + t231;
t248 = t110 * t109;
t191 = sqrt(-t248);
t198 = t137 * t166 - t138 * t167;
t277 = t198 * t191;
t179 = pkin(10) ^ 2;
t181 = pkin(8) ^ 2;
t187 = pkin(3) ^ 2;
t228 = t186 + t187;
t221 = -t181 + t228;
t115 = t121 + t179 + t221;
t120 = pkin(3) * t123 + pkin(4);
t100 = pkin(3) * t115 * t198 + t120 * t191;
t170 = cos(pkin(17));
t118 = t121 + t228;
t117 = 0.1e1 / t118 ^ 2;
t260 = pkin(4) * t117;
t226 = pkin(3) * t260;
t211 = t170 * t226;
t169 = sin(pkin(17));
t212 = t169 * t226;
t116 = 0.1e1 / t118;
t180 = 0.1e1 / pkin(10);
t246 = t116 * t180;
t266 = t169 / 0.2e1;
t98 = -pkin(3) * t277 + t115 * t120;
t91 = (-t98 * t170 / 0.2e1 + t100 * t266) * t246;
t90 = 0.1e1 / t91 ^ 2;
t92 = (t100 * t170 / 0.2e1 + t98 * t266) * t246;
t255 = t90 * t92;
t89 = 0.1e1 / t91;
t276 = (t100 * t211 + t98 * t212) * t89 - (t100 * t212 - t98 * t211) * t255;
t275 = -0.2e1 * t164 ^ 2;
t274 = pkin(3) * pkin(4);
t269 = t116 / 0.2e1;
t233 = t163 * t192;
t227 = pkin(1) * t256;
t239 = 0.1e1 / t192 * (t150 + t151) * t227;
t129 = (t233 + (-0.2e1 * t159 * pkin(1) - t155 - t239) * t164) * pkin(5);
t268 = -t129 / 0.2e1;
t131 = t159 * t239 + t185 * pkin(1) * t275 + (-t155 * t163 - t232) * pkin(5);
t267 = t131 / 0.2e1;
t171 = sin(pkin(19));
t265 = t171 / 0.2e1;
t264 = t173 / 0.2e1;
t178 = cos(pkin(15));
t263 = t178 / 0.2e1;
t259 = pkin(4) * t198;
t254 = t180 / (t90 * t92 ^ 2 + 0.1e1);
t213 = 0.1e1 / t158 ^ 2 * t227;
t107 = ((t129 * t264 + t177 * t267) * t156 + (t240 + t243) * t213) * t189;
t108 = ((t131 * t264 + t177 * t268) * t156 + (t241 - t242) * t213) * t189;
t103 = -t107 * t166 - t108 * t167;
t253 = t103 * t191;
t106 = 0.1e1 / t191;
t113 = t179 - t221 - 0.2e1 * t225;
t252 = t106 / t113;
t119 = -pkin(3) - t258;
t251 = t106 * t119;
t250 = t106 * t120;
t249 = t106 * t198;
t247 = t116 * t170;
t172 = cos(pkin(19));
t237 = t156 * t172;
t176 = sin(pkin(15));
t236 = t156 * t176;
t184 = 0.1e1 / pkin(6);
t235 = t156 * t184;
t222 = -t188 + t229;
t220 = t116 * t266;
t219 = -t247 / 0.2e1;
t218 = t247 / 0.2e1;
t217 = t156 * t265;
t216 = t156 * t263;
t215 = -0.2e1 * t187 * t259;
t214 = -0.2e1 * pkin(4) * t120 - t115;
t114 = -t179 + t181 + t118;
t97 = -pkin(4) * t277 - t114 * t119;
t96 = 0.1e1 / t97 ^ 2;
t99 = t114 * t259 - t119 * t191;
t210 = t118 / (t96 * t99 ^ 2 + 0.1e1);
t209 = (t109 + t110) * t274;
t112 = 0.1e1 / t113 ^ 2;
t208 = 0.2e1 * t112 * t191 * t274;
t205 = 0.1e1 / t97 * t210;
t154 = t161 + t183 + t222;
t160 = pkin(1) * t163 - pkin(5);
t142 = -pkin(1) * t232 - t154 * t160;
t202 = t142 * t213;
t201 = t143 * t213;
t200 = t144 * t213;
t145 = pkin(1) * t164 * t154 - t160 * t192;
t199 = t145 * t213;
t197 = pkin(4) * t96 * t99 * t210;
t102 = -t107 * t167 + t108 * t166;
t94 = t103 * t209;
t196 = -t102 * t191 - t249 * t94;
t101 = t198 * t209;
t195 = -t101 * t249 - t123 * t191;
t194 = pkin(3) * (t116 * t119 + t117 * t97);
t193 = pkin(3) * (-t116 * t186 * t198 + t260 * t99);
t153 = t183 - t222 + 0.2e1 * t224;
t152 = 0.1e1 / t153 ^ 2;
t139 = (t145 * t263 - t142 * t176 / 0.2e1) * t235;
t136 = (t142 * t263 + t145 * t176 / 0.2e1) * t235;
t135 = (t144 * t172 / 0.2e1 + t143 * t265) * t234;
t134 = (-t143 * t172 / 0.2e1 + t144 * t265) * t234;
t133 = 0.1e1 / t136 ^ 2;
t132 = 0.1e1 / t134 ^ 2;
t130 = -t160 * t239 + t190 * pkin(5) * t275 + (-t154 * t163 - t232) * pkin(1);
t128 = (t233 + (0.2e1 * t160 * pkin(5) - t154 - t239) * t164) * pkin(1);
t104 = 0.1e1 / (-t112 * t248 + 0.1e1);
t88 = t101 * t250 + t198 * t215 + (t115 * t123 - t277) * pkin(3);
t87 = (t198 * t214 + t195) * pkin(3);
t85 = t94 * t250 + t103 * t215 + (t102 * t115 - t253) * pkin(3);
t84 = (t103 * t214 + t196) * pkin(3);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, ((t85 * t218 + t84 * t220) * t89 - (t84 * t219 + t85 * t220) * t255 + t276 * t103) * t254, ((t88 * t218 + t87 * t220) * t89 - (t87 * t219 + t88 * t220) * t255 + t276 * t198) * t254, 0; 0, 0, 0, 1; 0, ((t130 * t216 + t178 * t199 - t128 * t236 / 0.2e1 - t176 * t202) / t136 - (t128 * t216 + t178 * t202 + t130 * t236 / 0.2e1 + t176 * t199) * t139 * t133) / (t133 * t139 ^ 2 + 0.1e1) * t184, 0, 0; 0, ((t129 * t217 + t171 * t201 + t172 * t200 + t237 * t267) / t134 - (t131 * t217 + t171 * t200 - t172 * t201 + t237 * t268) * t135 * t132) / (t132 * t135 ^ 2 + 0.1e1) * t189, 0, 0; 0, 0.2e1 * ((-t94 * t251 + (t102 * t114 - t253) * pkin(4)) * t269 + t103 * t193) * t205 - 0.2e1 * ((-t103 * t114 + t196) * t269 + t103 * t194) * t197, 0.2e1 * ((-t101 * t251 + (t123 * t114 - t277) * pkin(4)) * t269 + t198 * t193) * t205 - 0.2e1 * ((-t114 * t198 + t195) * t269 + t198 * t194) * t197, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, (-0.1e1 / t153 * t239 + 0.2e1 * pkin(1) * t152 * t223) / (-t152 * t238 + 0.1e1), 0, 0; 0, (t103 * t208 - t252 * t94) * t104, (-t101 * t252 + t198 * t208) * t104, 0;];
W = t1;
