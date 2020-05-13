% Calculate Gravitation load on the joints for
% palh3m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m1DE1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:30:11
% EndTime: 2020-04-18 10:41:49
% DurationCPUTime: 116.29s
% Computational Cost: add. (3241903->319), mult. (4894402->530), div. (216516->22), fcn. (3103493->36), ass. (0->255)
t317 = -m(8) - m(4);
t330 = m(9) - t317;
t322 = (pkin(3) * pkin(4));
t342 = 2 * t322;
t341 = pkin(3) * m(9);
t184 = pkin(4) ^ 2;
t183 = pkin(5) ^ 2;
t188 = pkin(1) ^ 2;
t171 = sin(qJ(2));
t173 = sin(pkin(16));
t177 = cos(qJ(2));
t309 = cos(pkin(16));
t149 = t171 * t173 - t177 * t309;
t307 = pkin(5) * t149;
t325 = -2 * pkin(1);
t277 = t307 * t325 + t188;
t136 = t183 + t277;
t273 = pkin(2) ^ 2 - pkin(6) ^ 2;
t125 = t136 + t273;
t141 = pkin(1) - t307;
t320 = -pkin(6) - pkin(2);
t122 = (pkin(5) - t320) * (pkin(5) + t320) + t277;
t321 = pkin(2) - pkin(6);
t123 = (pkin(5) - t321) * (pkin(5) + t321) + t277;
t189 = sqrt(-t123 * t122);
t151 = t171 * t309 + t177 * t173;
t306 = pkin(5) * t151;
t117 = t125 * t306 + t141 * t189;
t308 = sin(qJ(3));
t253 = t117 * t308;
t130 = 0.1e1 / t136;
t187 = 0.1e1 / pkin(2);
t285 = t130 * t187;
t283 = t151 * t189;
t116 = -pkin(5) * t283 + t125 * t141;
t176 = cos(qJ(3));
t291 = t116 * t176;
t110 = (-t291 / 0.2e1 + t253 / 0.2e1) * t285;
t254 = t116 * t308;
t290 = t117 * t176;
t111 = (t290 / 0.2e1 + t254 / 0.2e1) * t285;
t163 = pkin(18) + pkin(19);
t159 = sin(t163);
t160 = cos(t163);
t89 = t110 * t160 + t111 * t159;
t216 = pkin(4) * t89;
t292 = 0.2e1 * pkin(3) * t216 + t184;
t319 = -pkin(8) - pkin(10);
t80 = (pkin(3) - t319) * (pkin(3) + t319) + t292;
t318 = pkin(10) - pkin(8);
t81 = (pkin(3) - t318) * (pkin(3) + t318) + t292;
t190 = sqrt(-t81 * t80);
t71 = 0.1e1 / t190;
t340 = -t71 / 0.2e1;
t339 = -m(6) - m(5);
t165 = cos(pkin(17));
t185 = pkin(3) ^ 2;
t85 = t185 + t292;
t271 = 0.1e1 / t85 ^ 2 * t322;
t245 = t165 * t271;
t270 = pkin(1) * t306;
t289 = 0.2e1 / t189 * (t122 + t123) * t270;
t251 = -t289 / 0.2e1;
t284 = t149 * t189;
t102 = (t284 + (t141 * t325 - t125 + t251) * t151) * pkin(5);
t244 = 0.1e1 / t136 ^ 2 * t270;
t266 = t308 / 0.2e1;
t323 = -0.2e1 * t151 ^ 2;
t104 = t141 * t289 / 0.2e1 + t183 * pkin(1) * t323 + (-t125 * t149 - t283) * pkin(5);
t313 = t104 / 0.2e1;
t78 = ((t102 * t266 + t176 * t313) * t130 + (t254 + t290) * t244) * t187;
t314 = -t102 / 0.2e1;
t79 = ((t104 * t266 + t176 * t314) * t130 + (t253 - t291) * t244) * t187;
t69 = -t159 * t78 - t160 * t79;
t239 = t69 * t245;
t164 = sin(pkin(17));
t246 = t164 * t271;
t241 = t69 * t246;
t83 = 0.1e1 / t85;
t299 = t165 * t83;
t257 = t299 / 0.2e1;
t258 = -t299 / 0.2e1;
t312 = t164 / 0.2e1;
t259 = t83 * t312;
t180 = 0.1e1 / pkin(10);
t294 = t180 * t83;
t225 = t110 * t159 - t111 * t160;
t293 = t190 * t225;
t272 = -pkin(8) ^ 2 + pkin(10) ^ 2;
t82 = t85 + t272;
t86 = pkin(3) * t89 + pkin(4);
t65 = -pkin(3) * t293 + t82 * t86;
t315 = pkin(3) * t225;
t66 = t190 * t86 + t315 * t82;
t59 = (-t65 * t165 / 0.2e1 + t66 * t312) * t294;
t58 = 0.1e1 / t59 ^ 2;
t60 = (t66 * t165 / 0.2e1 + t65 * t312) * t294;
t295 = t180 / (t58 * t60 ^ 2 + 0.1e1);
t303 = t58 * t60;
t262 = t225 * t340;
t224 = (t80 + t81) * t342;
t62 = t69 * t224;
t68 = t159 * t79 - t160 * t78;
t219 = -t68 * t190 + t262 * t62;
t324 = -0.2e1 * pkin(4);
t255 = t324 * t86 - t82;
t31 = (t255 * t69 + t219) * pkin(3);
t248 = t185 * t225 * t324;
t263 = t71 * t86 / 0.2e1;
t32 = t62 * t263 + t69 * t248 + (-t190 * t69 + t68 * t82) * pkin(3);
t57 = 0.1e1 / t59;
t15 = ((t239 * t66 + t241 * t65 + t257 * t32 + t259 * t31) * t57 - (-t239 * t65 + t241 * t66 + t258 * t31 + t259 * t32) * t303) * t295;
t338 = t15 + 0.1e1;
t238 = t225 * t245;
t240 = t225 * t246;
t67 = t225 * t224;
t218 = -t89 * t190 + t262 * t67;
t54 = (t225 * t255 + t218) * pkin(3);
t55 = t67 * t263 + t225 * t248 + (t82 * t89 - t293) * pkin(3);
t18 = ((t238 * t66 + t240 * t65 + t257 * t55 + t259 * t54) * t57 - (-t238 * t65 + t240 * t66 + t258 * t54 + t259 * t55) * t303) * t295;
t337 = t18 + 0.1e1;
t302 = t83 * t85;
t170 = sin(qJ(4));
t175 = cos(qJ(4));
t336 = -m(6) * pkin(9) - t175 * mrSges(6,1) + t170 * mrSges(6,2) - mrSges(5,1);
t172 = sin(qJ(1));
t217 = t171 * t176 + t177 * t308;
t137 = t217 * t172;
t252 = t171 * t308;
t148 = -t177 * t176 + t252;
t138 = t148 * t172;
t47 = atan2(t60, t59);
t45 = cos(t47);
t35 = t138 * t45;
t44 = sin(t47);
t335 = -t137 * t44 - t35;
t168 = cos(pkin(19));
t166 = sin(pkin(19));
t311 = t166 / 0.2e1;
t107 = (-t116 * t168 / 0.2e1 + t117 * t311) * t285;
t108 = (t117 * t168 / 0.2e1 + t116 * t311) * t285;
t99 = atan2(t108, t107);
t93 = sin(t99);
t296 = t171 * t93;
t94 = cos(t99);
t92 = t177 * t94;
t76 = t92 - t296;
t214 = -t216 - pkin(3);
t230 = t85 - t272;
t209 = -pkin(4) * t293 - t214 * t230;
t208 = 0.1e1 / t209 ^ 2;
t220 = pkin(4) * t230;
t210 = -t190 * t214 + t220 * t225;
t207 = 0.1e1 / (t208 * t210 ^ 2 + 0.1e1);
t201 = t207 * t208 * t210;
t202 = t207 / t209;
t211 = t214 * t342 - t220;
t213 = t214 * t340;
t223 = -pkin(4) * t190 - 0.2e1 * t184 * t315;
t334 = (t213 * t67 + t220 * t89 + t223 * t225) * t202 - (pkin(4) * t218 + t211 * t225) * t201;
t333 = (t213 * t62 + t220 * t68 + t223 * t69) * t202 - (pkin(4) * t219 + t211 * t69) * t201;
t167 = sin(pkin(18));
t169 = cos(pkin(18));
t267 = t83 / pkin(8) / 0.2e1;
t206 = atan2(t210 * t267, t209 * t267);
t205 = sin(t206);
t61 = cos(t206);
t332 = -t167 * t61 + t169 * t205;
t43 = -t167 * t205 - t169 * t61;
t329 = t302 * t43;
t243 = -m(6) * pkin(11) + mrSges(5,2) - mrSges(6,3);
t328 = t332 * t302;
t327 = t336 * t338;
t326 = t336 * t337;
t105 = 0.1e1 / t107 ^ 2;
t232 = t116 * t244;
t235 = t117 * t244;
t250 = t130 * t311;
t288 = t130 * t168;
t63 = ((t102 * t250 + t166 * t232 + t168 * t235 + t288 * t313) / t107 - (t104 * t250 + t166 * t235 - t168 * t232 + t288 * t314) * t108 * t105) / (t105 * t108 ^ 2 + 0.1e1) * t187;
t316 = t63 + 0.1e1;
t179 = cos(pkin(15));
t310 = t179 / 0.2e1;
t132 = t137 * pkin(4);
t178 = cos(qJ(1));
t140 = t217 * t178;
t134 = t140 * pkin(4);
t146 = t148 * pkin(4);
t162 = t177 * pkin(1);
t280 = t177 * t178;
t139 = t176 * t280 - t178 * t252;
t227 = t139 * t44 + t140 * t45;
t226 = t148 * t45 + t217 * t44;
t301 = mrSges(6,2) * t175;
t33 = t137 * t45;
t37 = t139 * t45;
t41 = t217 * t45;
t174 = sin(pkin(15));
t287 = t130 * t174;
t182 = 0.1e1 / pkin(6);
t286 = t130 * t182;
t282 = t171 * t172;
t279 = t137 * mrSges(4,1) - t138 * mrSges(4,2);
t278 = t140 * mrSges(4,1) + t139 * mrSges(4,2);
t276 = t148 * mrSges(4,1) + mrSges(4,2) * t217;
t274 = pkin(1) * t280 + t178 * pkin(13);
t269 = pkin(1) * t282;
t268 = pkin(1) * t171 * t178;
t265 = t338 * t44;
t264 = t337 * t44;
t261 = t171 * t316;
t260 = t177 * t316;
t256 = -t139 * pkin(4) + t274;
t249 = t130 * t310;
t242 = t94 * t260;
t234 = t174 * t244;
t233 = t179 * t244;
t124 = t136 - t273;
t142 = pkin(1) * t149 - pkin(5);
t115 = -pkin(1) * t283 - t124 * t142;
t118 = pkin(1) * t151 * t124 - t142 * t189;
t109 = (t115 * t310 + t118 * t174 / 0.2e1) * t286;
t112 = (t118 * t310 - t115 * t174 / 0.2e1) * t286;
t100 = atan2(t112, t109);
t96 = sin(t100);
t97 = cos(t100);
t231 = mrSges(7,1) * t97 - mrSges(7,2) * t96;
t229 = mrSges(3,1) * t177 - mrSges(3,2) * t171;
t77 = t171 * t94 + t177 * t93;
t221 = -mrSges(8,3) + mrSges(2,2) - mrSges(7,3) - mrSges(3,3) - mrSges(9,3) - mrSges(5,3) - mrSges(4,3);
t52 = -t260 * t93 - t261 * t94;
t212 = pkin(7) * m(7) - m(3) * pkin(13) - mrSges(2,1) - t229 - t231;
t158 = t178 * t170;
t106 = 0.1e1 / t109 ^ 2;
t103 = t142 * t251 + t188 * pkin(5) * t323 + (-t124 * t149 - t283) * pkin(1);
t101 = (t284 + (0.2e1 * t142 * pkin(5) - t124 + t251) * t151) * pkin(1);
t91 = t93 * t282;
t75 = t76 * t178;
t74 = t77 * t178;
t73 = t172 * t92 - t91;
t72 = t77 * t172;
t64 = ((t103 * t249 + t118 * t233 - t101 * t287 / 0.2e1 - t115 * t234) / t109 - (t101 * t249 + t115 * t233 + t103 * t287 / 0.2e1 + t118 * t234) * t112 * t106) / (t106 * t112 ^ 2 + 0.1e1) * t182;
t53 = -t261 * t93 + t63 * t92 + t92;
t51 = t52 * t178;
t50 = (t296 * t316 - t242) * t178;
t49 = t52 * t172;
t48 = t91 + (t296 * t63 - t242) * t172;
t29 = t140 * t44 - t37;
t24 = t170 * t172 + t175 * t29;
t23 = -t170 * t29 + t172 * t175;
t17 = t334 * t328;
t16 = t334 * t329;
t14 = t333 * t328;
t13 = t333 * t329;
t1 = [(-(t332 * t74 + t43 * t75) * mrSges(9,1) - (t332 * t75 - t43 * t74) * mrSges(9,2) - mrSges(8,1) * t75 + mrSges(8,2) * t74 - m(9) * ((t167 * t74 + t169 * t75) * pkin(3) + t274) - t24 * mrSges(6,1) + mrSges(4,1) * t139 - mrSges(4,2) * t140 - m(6) * (pkin(9) * t29 + t256) - m(5) * t256 - t29 * mrSges(5,1) - t23 * mrSges(6,2) - t243 * t227 + t317 * t274 + t212 * t178 + t221 * t172) * g(2) + (-(-t332 * t72 - t43 * t73) * mrSges(9,1) - (-t332 * t73 + t43 * t72) * mrSges(9,2) + mrSges(8,1) * t73 - mrSges(8,2) * t72 - t158 * mrSges(6,1) + mrSges(4,2) * t137 + t336 * t335 + t243 * t33 - (-t167 * t72 - t169 * t73) * t341 + (t221 - t301) * t178 + (-t212 + (-t339 + t330) * (pkin(13) + t162)) * t172 - (t339 * pkin(4) + t243 * t44 - mrSges(4,1)) * t138) * g(1), (-(-t167 * t48 + t169 * t49) * t341 - (-t13 * t72 + t14 * t73 - t332 * t48 + t43 * t49) * mrSges(9,1) - (-t13 * t73 - t14 * t72 + t332 * t49 + t43 * t48) * mrSges(9,2) - t279 - mrSges(8,1) * t49 - mrSges(8,2) * t48 + t330 * t269 + t339 * (t132 - t269) + t336 * (-t138 * t265 + t15 * t33 + t33) + t243 * (t137 * t265 + t15 * t35 + t35)) * g(2) + (-(-t167 * t50 + t169 * t51) * t341 - (-t13 * t74 + t14 * t75 - t332 * t50 + t43 * t51) * mrSges(9,1) - (-t13 * t75 - t14 * t74 + t332 * t51 + t43 * t50) * mrSges(9,2) - t278 - mrSges(8,1) * t51 - mrSges(8,2) * t50 + t227 * t327 + t339 * (t134 - t268) + t243 * (t140 * t265 - t15 * t37 - t37) + t330 * t268) * g(1) + (-(-t167 * t52 + t169 * t53) * t341 - (t13 * t76 + t14 * t77 - t332 * t52 + t43 * t53) * mrSges(9,1) - (-t13 * t77 + t14 * t76 + t332 * t53 + t43 * t52) * mrSges(9,2) - t276 - mrSges(8,1) * t53 - mrSges(8,2) * t52 - t231 * t64 - t229 + t226 * t327 + t339 * (t146 + t162) + t243 * (t148 * t265 - t15 * t41 - t41) - t330 * t162) * g(3) + (mrSges(3,1) * t171 + mrSges(3,2) * t177 - (-mrSges(7,1) * t96 - mrSges(7,2) * t97) * t64) * (g(1) * t178 + g(2) * t172), (-(t16 * t76 + t17 * t77) * mrSges(9,1) - (-t16 * t77 + t17 * t76) * mrSges(9,2) - t276 + t339 * t146 + t226 * t326 + t243 * (t148 * t264 - t18 * t41 - t41)) * g(3) + (-(-t16 * t72 + t17 * t73) * mrSges(9,1) - (-t16 * t73 - t17 * t72) * mrSges(9,2) - t279 + t336 * (-t138 * t264 + t18 * t33 + t33) + t339 * t132 + t243 * (t137 * t264 + t18 * t35 + t35)) * g(2) + (-(-t16 * t74 + t17 * t75) * mrSges(9,1) - (-t16 * t75 - t17 * t74) * mrSges(9,2) - t278 + t339 * t134 + t243 * (t140 * t264 - t18 * t37 - t37) + t227 * t326) * g(1), -g(2) * ((t170 * t335 - t175 * t178) * mrSges(6,1) + (t175 * t335 + t158) * mrSges(6,2)) - g(1) * (mrSges(6,1) * t23 - mrSges(6,2) * t24) - g(3) * (-mrSges(6,1) * t170 - t301) * (t148 * t44 - t41)];
taug = t1(:);
