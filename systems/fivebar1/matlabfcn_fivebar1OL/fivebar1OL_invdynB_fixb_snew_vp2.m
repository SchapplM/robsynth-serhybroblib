% Calculate vector of inverse dynamics base forces with Newton-Euler for
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = fivebar1OL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1OL_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:12
% EndTime: 2020-04-27 06:13:12
% DurationCPUTime: 0.38s
% Computational Cost: add. (342->101), mult. (512->123), div. (0->0), fcn. (340->8), ass. (0->61)
t252 = sin(qJ(1));
t260 = pkin(2) * m(3);
t251 = sin(qJ(2));
t248 = t251 * mrSges(3,2);
t255 = cos(qJ(2));
t287 = t255 * mrSges(3,1);
t291 = -t287 + t248;
t304 = -mrSges(2,1) - t260 - t291;
t306 = t252 * t304;
t249 = sin(qJ(4));
t247 = t249 * mrSges(5,2);
t253 = cos(qJ(4));
t288 = t253 * mrSges(5,1);
t237 = -t247 + t288;
t290 = pkin(3) * m(5);
t228 = mrSges(4,1) + t237 + t290;
t236 = t249 * mrSges(5,1) + t253 * mrSges(5,2);
t231 = mrSges(4,2) + t236;
t250 = sin(qJ(3));
t254 = cos(qJ(3));
t276 = -t228 * t254 + t231 * t250;
t256 = cos(qJ(1));
t284 = t304 * t256;
t305 = -t276 + mrSges(1,1) - t284;
t262 = qJD(3) ^ 2;
t267 = t250 * t236 - t237 * t254;
t275 = t250 * t228 + t231 * t254;
t303 = -t276 * qJDD(3) - t267 * qJDD(4) - t275 * t262;
t239 = t251 * mrSges(3,1) + t255 * mrSges(3,2);
t297 = -mrSges(2,2) + t239;
t302 = t252 * t297;
t294 = (2 * qJD(3) + qJD(4)) * qJD(4);
t293 = (2 * qJD(1) + qJD(2)) * qJD(2);
t225 = t297 * t256;
t300 = -mrSges(1,2) + t225 - t275;
t292 = t284 - t302;
t259 = m(4) + m(5);
t235 = t253 * Ifges(5,5) - t249 * Ifges(5,6);
t274 = t225 + t306;
t238 = t255 * Ifges(3,5) - t251 * Ifges(3,6);
t241 = Ifges(3,5) * t251 + Ifges(3,6) * t255;
t234 = t249 * Ifges(5,5) + t253 * Ifges(5,6);
t226 = -pkin(3) * mrSges(5,3) + Ifges(4,5) + t235;
t230 = Ifges(4,6) + t234;
t271 = t250 * t226 + t230 * t254;
t227 = -pkin(2) * mrSges(3,3) + Ifges(2,5) - t238;
t233 = -Ifges(2,6) + t241;
t270 = t252 * t227 - t233 * t256;
t269 = t234 * t254 + t250 * t235;
t214 = -t250 * t234 + t235 * t254;
t268 = -t254 * t236 - t250 * t237;
t215 = -t238 * t256 + t252 * t241;
t266 = t252 * t238 + t241 * t256;
t216 = t239 * t256 - t252 * t291;
t213 = t252 * t239 + t256 * t291;
t264 = qJD(1) ^ 2;
t245 = m(2) + m(3) + m(1) + t259;
t244 = mrSges(2,3) + mrSges(4,3) + mrSges(5,3) + mrSges(3,3) + mrSges(1,3);
t212 = t227 * t256 + t252 * t233;
t211 = t226 * t254 - t250 * t230;
t1 = [-t245 * g(1) + t274 * qJDD(1) + t216 * qJDD(2) - t275 * qJDD(3) + t268 * qJDD(4) - t293 * t213 + t276 * t262 + t292 * t264 + t294 * t267; -t245 * g(2) - t292 * qJDD(1) + t213 * qJDD(2) + t293 * t216 + t274 * t264 + t294 * t268 + t303; -t245 * g(3); -t271 * t262 + t211 * qJDD(3) + t214 * qJDD(4) - (-t300 - t306) * g(3) + t215 * qJDD(2) - t270 * t264 + t212 * qJDD(1) + t244 * g(2) - t294 * t269 + t293 * t266; t269 * qJDD(4) + t211 * t262 + t271 * qJDD(3) - (-pkin(1) * t259 - t302 - t305) * g(3) - t244 * g(1) - t266 * qJDD(2) + t212 * t264 + t270 * qJDD(1) + t293 * t215 + t294 * t214; Ifges(5,3) * qJDD(4) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + Ifges(3,3) * qJDD(2) - t300 * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1) - t305 * g(2) + (-g(1) * t304 - g(2) * t297) * t252 + t294 * (t268 * pkin(1) - pkin(3) * t236) + (t237 * qJDD(4) + (-0.2e1 * t247 + 0.2e1 * t288 + t290) * qJDD(3)) * pkin(3) + (-t259 * g(2) + t303) * pkin(1) + (t291 * qJDD(2) + t293 * t239 + (0.2e1 * t248 - 0.2e1 * t287 + t260) * qJDD(1)) * pkin(2);];
tauB = t1;
