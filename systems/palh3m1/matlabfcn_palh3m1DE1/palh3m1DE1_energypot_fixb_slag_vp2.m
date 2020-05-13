% Calculate potential energy for
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1DE1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_energypot_fixb_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_energypot_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE1_energypot_fixb_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:25:40
% EndTime: 2020-04-18 10:25:50
% DurationCPUTime: 7.09s
% Computational Cost: add. (157249->148), mult. (236353->206), div. (11232->6), fcn. (149713->36), ass. (0->104)
t238 = cos(qJ(2));
t294 = pkin(1) * t238 + pkin(13);
t293 = pkin(3) * m(9);
t292 = -m(6) - m(5);
t291 = -m(4) - m(8) - m(9);
t232 = sin(qJ(2));
t234 = sin(pkin(16));
t240 = cos(pkin(16));
t210 = t232 * t234 - t238 * t240;
t274 = pkin(5) * t210;
t267 = (-0.2e1 * t274 + pkin(1)) * pkin(1);
t260 = pkin(5) ^ 2 + t267;
t265 = pkin(2) ^ 2 - pkin(6) ^ 2;
t196 = t260 - t265;
t206 = pkin(1) * t210 - pkin(5);
t285 = -pkin(6) + pkin(2);
t286 = -pkin(6) - pkin(2);
t195 = sqrt(-((pkin(5) - t285) * (pkin(5) + t285) + t267) * ((pkin(5) - t286) * (pkin(5) + t286) + t267));
t212 = t232 * t240 + t234 * t238;
t271 = t195 * t212;
t191 = -pkin(1) * t271 - t196 * t206;
t194 = pkin(1) * t196 * t212 - t195 * t206;
t235 = sin(pkin(15));
t198 = 0.1e1 / t260;
t270 = t198 / pkin(6);
t277 = cos(pkin(15)) / 0.2e1;
t188 = atan2((t194 * t277 - t191 * t235 / 0.2e1) * t270, (t191 * t277 + t194 * t235 / 0.2e1) * t270);
t185 = sin(t188);
t186 = cos(t188);
t290 = pkin(7) * m(7) - m(3) * pkin(13) - mrSges(3,1) * t238 - mrSges(7,1) * t186 + mrSges(3,2) * t232 + mrSges(7,2) * t185 - mrSges(2,1);
t289 = -m(6) * pkin(11) + mrSges(5,2) - mrSges(6,3);
t230 = sin(qJ(4));
t236 = cos(qJ(4));
t288 = -m(6) * pkin(9) - t236 * mrSges(6,1) + t230 * mrSges(6,2) - mrSges(5,1);
t287 = t230 * mrSges(6,1) + t236 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3);
t284 = -pkin(8) - pkin(10);
t283 = -pkin(8) + pkin(10);
t197 = t260 + t265;
t205 = pkin(1) - t274;
t192 = -pkin(5) * t271 + t197 * t205;
t282 = -t192 / 0.2e1;
t193 = pkin(5) * t197 * t212 + t195 * t205;
t281 = t193 / 0.2e1;
t280 = sin(pkin(17)) / 0.2e1;
t279 = sin(pkin(19)) / 0.2e1;
t231 = sin(qJ(3));
t278 = t231 / 0.2e1;
t237 = cos(qJ(3));
t269 = t198 / pkin(2);
t189 = (t193 * t278 + t237 * t282) * t269;
t190 = (t192 * t278 + t237 * t281) * t269;
t223 = pkin(18) + pkin(19);
t218 = sin(t223);
t219 = cos(t223);
t182 = -t189 * t219 - t190 * t218;
t275 = pkin(4) * t182;
t217 = t232 * pkin(1) + pkin(12);
t268 = -0.2e1 * pkin(3) * t275 + pkin(4) ^ 2;
t168 = sqrt(-((pkin(3) - t283) * (pkin(3) + t283) + t268) * ((pkin(3) - t284) * (pkin(3) + t284) + t268));
t181 = t189 * t218 - t190 * t219;
t273 = t168 * t181;
t263 = pkin(3) ^ 2 + t268;
t177 = 0.1e1 / t263;
t272 = t177 / pkin(10);
t233 = sin(qJ(1));
t213 = t294 * t233;
t239 = cos(qJ(1));
t214 = t294 * t239;
t266 = pkin(8) ^ 2 - pkin(10) ^ 2;
t259 = t177 / pkin(8) / 0.2e1;
t228 = cos(pkin(19));
t187 = atan2((t192 * t279 + t228 * t281) * t269, (t193 * t279 + t228 * t282) * t269);
t183 = sin(t187);
t184 = cos(t187);
t174 = t183 * t238 + t184 * t232;
t173 = -t183 * t232 + t184 * t238;
t258 = t231 * t238 + t232 * t237;
t209 = t231 * t232 - t237 * t238;
t175 = t263 + t266;
t178 = -pkin(3) + t275;
t255 = atan2((pkin(4) * t175 * t181 - t168 * t178) * t259, (-pkin(4) * t273 - t175 * t178) * t259);
t254 = sin(t255);
t229 = cos(pkin(18));
t227 = sin(pkin(18));
t225 = cos(pkin(17));
t204 = t209 * t239;
t203 = t258 * t239;
t202 = t209 * t233;
t201 = t258 * t233;
t179 = -pkin(3) * t182 + pkin(4);
t176 = t263 - t266;
t172 = t173 * t239;
t171 = t174 * t239;
t170 = t173 * t233;
t169 = t174 * t233;
t167 = pkin(3) * t176 * t181 + t168 * t179;
t166 = -pkin(3) * t273 + t176 * t179;
t165 = cos(t255);
t163 = atan2((t167 * t225 / 0.2e1 + t166 * t280) * t272, (-t166 * t225 / 0.2e1 + t167 * t280) * t272);
t162 = cos(t163);
t161 = sin(t163);
t160 = -t229 * t165 - t227 * t254;
t159 = t165 * t227 - t229 * t254;
t1 = (-mrSges(8,1) * t174 - mrSges(8,2) * t173 + mrSges(4,1) * t258 - mrSges(4,2) * t209 - mrSges(1,3) - mrSges(2,3) - m(7) * pkin(14) - t185 * mrSges(7,1) - t186 * mrSges(7,2) - (-t173 * t227 + t174 * t229) * t293 - (t159 * t173 + t160 * t174) * mrSges(9,1) - (-t159 * t174 + t160 * t173) * mrSges(9,2) - mrSges(3,1) * t232 - mrSges(3,2) * t238 + t292 * (-pkin(4) * t258 + t217) + t291 * t217 + t288 * (t161 * t209 - t162 * t258) + t289 * (-t161 * t258 - t209 * t162) + (-m(2) - m(7) - m(3)) * pkin(12)) * g(3) + (-(-t159 * t169 + t160 * t170) * mrSges(9,1) - (-t159 * t170 - t160 * t169) * mrSges(9,2) - mrSges(8,1) * t170 + mrSges(8,2) * t169 - mrSges(4,1) * t202 - mrSges(4,2) * t201 - (t169 * t227 + t170 * t229) * t293 - mrSges(1,2) + t292 * (t202 * pkin(4) + t213) + t291 * t213 + t289 * (t161 * t202 - t201 * t162) + t290 * t233 + t288 * (t161 * t201 + t162 * t202) + t287 * t239) * g(2) + (-mrSges(4,1) * t204 - mrSges(4,2) * t203 - mrSges(8,1) * t172 + mrSges(8,2) * t171 - (-t159 * t171 + t160 * t172) * mrSges(9,1) - (-t159 * t172 - t160 * t171) * mrSges(9,2) - (t171 * t227 + t172 * t229) * t293 - mrSges(1,1) + t292 * (t204 * pkin(4) + t214) + t291 * t214 + t289 * (t161 * t204 - t203 * t162) + t290 * t239 + t288 * (t161 * t203 + t162 * t204) - t287 * t233) * g(1);
U = t1;
