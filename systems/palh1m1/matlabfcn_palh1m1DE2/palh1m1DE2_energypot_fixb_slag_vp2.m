% Calculate potential energy for
% palh1m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1DE2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_energypot_fixb_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE2_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE2_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:02:11
% EndTime: 2020-04-14 20:02:18
% DurationCPUTime: 4.38s
% Computational Cost: add. (70637->155), mult. (105949->176), div. (5052->9), fcn. (67021->48), ass. (0->118)
t326 = -m(5) - m(6);
t325 = pkin(4) * m(11);
t242 = sin(qJ(4));
t248 = cos(qJ(4));
t324 = -m(6) * pkin(10) - t248 * mrSges(6,1) + t242 * mrSges(6,2) - mrSges(5,1);
t323 = -m(9) - m(3) - m(10);
t322 = -m(4) - m(8) - m(11);
t321 = m(10) * pkin(2) + mrSges(9,1);
t320 = m(6) * pkin(12) - mrSges(5,2) + mrSges(6,3);
t319 = -mrSges(6,1) * t242 - mrSges(6,2) * t248 + mrSges(2,2) - mrSges(11,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3);
t270 = pkin(1) ^ 2;
t244 = sin(qJ(2));
t246 = sin(pkin(19));
t250 = cos(qJ(2));
t252 = cos(pkin(19));
t223 = t244 * t252 - t246 * t250;
t301 = pkin(7) * t223;
t293 = -0.2e1 * pkin(1) * t301 + t270;
t314 = -pkin(8) + pkin(3);
t315 = -pkin(8) - pkin(3);
t210 = sqrt(-((pkin(7) - t314) * (pkin(7) + t314) + t293) * ((pkin(7) - t315) * (pkin(7) + t315) + t293));
t286 = pkin(7) ^ 2 + t293;
t290 = pkin(3) ^ 2 - pkin(8) ^ 2;
t213 = t286 + t290;
t218 = pkin(1) - t301;
t224 = t244 * t246 + t250 * t252;
t207 = pkin(7) * t213 * t224 + t210 * t218;
t243 = sin(qJ(3));
t215 = 0.1e1 / t286;
t296 = t215 / pkin(3);
t249 = cos(qJ(3));
t307 = -t249 / 0.2e1;
t298 = t210 * t224;
t206 = -pkin(7) * t298 + t213 * t218;
t311 = -t206 / 0.2e1;
t203 = (t207 * t307 + t243 * t311) * t296;
t310 = t207 / 0.2e1;
t204 = (t206 * t307 + t243 * t310) * t296;
t234 = pkin(23) + pkin(22);
t229 = sin(t234);
t230 = cos(t234);
t185 = t203 * t230 + t204 * t229;
t304 = pkin(5) * t185;
t295 = -0.2e1 * pkin(4) * t304 + pkin(5) ^ 2;
t287 = pkin(4) ^ 2 + t295;
t292 = pkin(9) ^ 2 - pkin(11) ^ 2;
t179 = t287 - t292;
t182 = -pkin(4) * t185 + pkin(5);
t312 = -pkin(9) + pkin(11);
t313 = -pkin(9) - pkin(11);
t177 = sqrt(-((pkin(4) - t312) * (pkin(4) + t312) + t295) * ((pkin(4) - t313) * (pkin(4) + t313) + t295));
t186 = -t203 * t229 + t204 * t230;
t300 = t177 * t186;
t175 = -pkin(4) * t300 + t179 * t182;
t176 = pkin(4) * t179 * t186 + t177 * t182;
t235 = qJ(2) + qJ(3);
t240 = cos(pkin(21));
t180 = 0.1e1 / t287;
t299 = t180 / pkin(11);
t308 = sin(pkin(21)) / 0.2e1;
t171 = atan2((t175 * t308 + t176 * t240 / 0.2e1) * t299, (-t175 * t240 / 0.2e1 + t176 * t308) * t299) + t235;
t169 = sin(t171);
t170 = cos(t171);
t178 = t287 + t292;
t181 = -pkin(4) + t304;
t239 = cos(pkin(23));
t309 = sin(pkin(23)) / 0.2e1;
t190 = qJ(2) + atan2((t206 * t309 + t239 * t310) * t296, (t207 * t309 + t239 * t311) * t296);
t187 = pkin(22) - t190;
t284 = t180 / pkin(9) / 0.2e1;
t174 = -atan2((pkin(5) * t178 * t186 - t177 * t181) * t284, (-pkin(5) * t300 - t178 * t181) * t284) + t187;
t172 = sin(t174);
t173 = cos(t174);
t188 = sin(t190);
t189 = cos(t190);
t212 = t286 - t290;
t219 = pkin(1) * t223 - pkin(7);
t205 = -pkin(1) * t298 - t212 * t219;
t208 = pkin(1) * t212 * t224 - t210 * t219;
t247 = sin(pkin(18));
t297 = t215 / pkin(8);
t306 = cos(pkin(18)) / 0.2e1;
t194 = atan2((t208 * t306 + t205 * t247 / 0.2e1) * t297, (t205 * t306 - t247 * t208 / 0.2e1) * t297);
t191 = sin(t194);
t192 = cos(t194);
t238 = sin(pkin(20));
t241 = cos(pkin(20));
t303 = pkin(6) * (-t238 * t249 - t241 * t243);
t289 = pkin(1) * t303;
t217 = -0.2e1 * t289;
t263 = pkin(6) ^ 2;
t294 = t217 + t263;
t316 = -pkin(2) + pkin(13);
t317 = -pkin(2) - pkin(13);
t209 = sqrt(-((pkin(1) - t316) * (pkin(1) + t316) + t294) * ((pkin(1) - t317) * (pkin(1) + t317) + t294));
t268 = pkin(2) ^ 2;
t291 = t263 + t270;
t285 = -pkin(13) ^ 2 + t291;
t211 = t217 + t268 + t285;
t216 = -pkin(1) + t303;
t305 = 0.1e1 / pkin(2) / 0.2e1;
t283 = 0.1e1 / (t217 + t291) * t305;
t302 = pkin(6) * (t238 * t243 - t241 * t249);
t201 = qJ(2) + atan2((-t209 * t216 + t211 * t302) * t283, (-t209 * t302 - t211 * t216) * t283);
t282 = 0.1e1 / pkin(13) * t305;
t197 = atan2(t209 * t282, (t268 - t285 + 0.2e1 * t289) * t282) + t201;
t195 = sin(t197);
t196 = cos(t197);
t199 = sin(t201);
t200 = cos(t201);
t227 = -t244 * pkin(1) + pkin(16);
t231 = sin(t235);
t232 = cos(t235);
t318 = -sin(t187) * t325 + mrSges(9,2) * t200 + mrSges(3,1) * t244 + mrSges(3,2) * t250 - mrSges(7,1) * t192 + mrSges(7,2) * t191 + mrSges(11,1) * t172 - mrSges(11,2) * t173 + m(7) * pkin(15) - mrSges(10,1) * t195 - mrSges(10,2) * t196 + mrSges(8,1) * t188 + mrSges(8,2) * t189 - mrSges(4,1) * t232 + mrSges(4,2) * t231 - mrSges(2,1) + t322 * t227 + t326 * (pkin(5) * t232 + t227) + t321 * t199 + t324 * t170 - t320 * t169 + t323 * pkin(16);
t228 = t250 * pkin(1) + pkin(14);
t251 = cos(qJ(1));
t245 = sin(qJ(1));
t1 = (t173 * mrSges(11,1) - t195 * mrSges(10,2) + t196 * mrSges(10,1) + mrSges(9,2) * t199 + m(7) * pkin(17) - t191 * mrSges(7,1) - t192 * mrSges(7,2) - cos(t187) * t325 + t172 * mrSges(11,2) - mrSges(3,1) * t250 + mrSges(3,2) * t244 - mrSges(8,1) * t189 + mrSges(8,2) * t188 - mrSges(4,1) * t231 - mrSges(4,2) * t232 - mrSges(1,3) - mrSges(2,3) + t326 * (pkin(5) * t231 + t228) + t322 * t228 - t321 * t200 + t320 * t170 + t324 * t169 + (-m(2) - m(7) + t323) * pkin(14)) * g(3) + (t318 * t245 - t319 * t251 - mrSges(1,2)) * g(2) + (t319 * t245 + t318 * t251 - mrSges(1,1)) * g(1);
U = t1;
