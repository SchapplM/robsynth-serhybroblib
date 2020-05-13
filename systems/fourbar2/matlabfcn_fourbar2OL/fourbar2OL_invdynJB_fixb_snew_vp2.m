% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% fourbar2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = fourbar2OL_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2OL_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_invdynJB_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2OL_invdynJB_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar2OL_invdynJB_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar2OL_invdynJB_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:25
% EndTime: 2020-04-24 20:32:26
% DurationCPUTime: 0.25s
% Computational Cost: add. (236->83), mult. (355->103), div. (0->0), fcn. (220->6), ass. (0->46)
t217 = sin(qJ(1));
t222 = m(3) * pkin(2);
t216 = sin(qJ(2));
t213 = t216 * mrSges(3,2);
t219 = cos(qJ(2));
t240 = t219 * mrSges(3,1);
t244 = -t240 + t213;
t252 = -mrSges(2,1) - t222 - t244;
t253 = t217 * t252;
t205 = t216 * mrSges(3,1) + t219 * mrSges(3,2);
t248 = -mrSges(2,2) + t205;
t251 = t217 * t248;
t246 = (2 * qJD(1) + qJD(2)) * qJD(2);
t215 = sin(qJ(3));
t218 = cos(qJ(3));
t202 = t215 * mrSges(4,1) + mrSges(4,2) * t218;
t208 = mrSges(4,1) * t218 - t215 * mrSges(4,2);
t223 = qJD(3) ^ 2;
t250 = t208 * qJDD(3) - t202 * t223;
t220 = cos(qJ(1));
t193 = t248 * t220;
t247 = t193 - mrSges(1,2) - t202;
t238 = t252 * t220;
t245 = -t238 + mrSges(1,1) + t208;
t183 = -t238 + t251;
t235 = Ifges(4,3) * qJDD(3);
t184 = t193 + t253;
t230 = (Ifges(2,3) + Ifges(3,3) + (0.2e1 * t213 - 0.2e1 * t240 + t222) * pkin(2)) * qJDD(1) + (t244 * pkin(2) + Ifges(3,3)) * qJDD(2) + t246 * pkin(2) * t205;
t204 = t219 * Ifges(3,5) - t216 * Ifges(3,6);
t207 = Ifges(3,5) * t216 + Ifges(3,6) * t219;
t194 = -mrSges(3,3) * pkin(2) + Ifges(2,5) - t204;
t200 = -Ifges(2,6) + t207;
t228 = t217 * t194 - t200 * t220;
t186 = -t204 * t220 + t217 * t207;
t227 = t217 * t204 + t207 * t220;
t187 = t205 * t220 - t217 * t244;
t185 = t217 * t205 + t220 * t244;
t225 = qJD(1) ^ 2;
t214 = m(2) + m(3) + m(4) + m(1);
t211 = mrSges(1,3) + mrSges(2,3) + mrSges(3,3) + mrSges(4,3);
t203 = t218 * Ifges(4,5) - t215 * Ifges(4,6);
t201 = -t215 * Ifges(4,5) - t218 * Ifges(4,6);
t198 = -t220 * g(1) - t217 * g(2) - t225 * pkin(2);
t197 = t217 * g(1) - t220 * g(2) + qJDD(1) * pkin(2);
t182 = t194 * t220 + t217 * t200;
t1 = [-t214 * g(1) + t184 * qJDD(1) + t187 * qJDD(2) - t202 * qJDD(3) - t183 * t225 - t246 * t185 - t208 * t223; -t214 * g(2) + t183 * qJDD(1) + t185 * qJDD(2) + t184 * t225 + t246 * t187 + t250; -t214 * g(3); t211 * g(2) - (-t247 - t253) * g(3) + t182 * qJDD(1) - t228 * t225 + t186 * qJDD(2) + t203 * qJDD(3) + t201 * t223 + t246 * t227; -t211 * g(1) - (-m(4) * pkin(1) - t245 - t251) * g(3) + t228 * qJDD(1) + t182 * t225 - t227 * qJDD(2) - t201 * qJDD(3) + t203 * t223 + t246 * t186; -t245 * g(2) - t247 * g(1) + t235 + (-g(1) * t252 - g(2) * t248) * t217 + (-m(4) * g(2) + t250) * pkin(1) + t230; -t184 * g(1) - t183 * g(2) + t230; Ifges(3,3) * (qJDD(1) + qJDD(2)) + mrSges(3,1) * (-t219 * t197 + t216 * t198) - mrSges(3,2) * (-t216 * t197 - t219 * t198); t202 * g(1) - t208 * g(2) + t235; 0;];
tauJB = t1;
