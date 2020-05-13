% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% fourbarprisOL
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
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
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = fourbarprisOL_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisOL_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_invdynJB_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_invdynJB_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisOL_invdynJB_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisOL_invdynJB_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:06
% EndTime: 2020-05-07 09:52:06
% DurationCPUTime: 0.11s
% Computational Cost: add. (207->75), mult. (323->87), div. (0->0), fcn. (104->4), ass. (0->30)
t142 = -m(2) - m(3);
t141 = -mrSges(2,1) - mrSges(3,3);
t140 = Ifges(3,4) + Ifges(2,6);
t139 = Ifges(2,5) - Ifges(3,6);
t129 = sin(qJ(1));
t131 = cos(qJ(1));
t119 = t131 * g(1) + t129 * g(2);
t133 = qJD(1) ^ 2;
t113 = -t133 * qJ(2) + qJDD(2) + t119;
t138 = m(3) * t113 + qJDD(1) * mrSges(3,1);
t117 = -t129 * g(1) + t131 * g(2);
t110 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t117;
t137 = -m(3) * t110 + t133 * mrSges(3,1) + qJDD(1) * mrSges(3,3);
t128 = sin(qJ(3));
t130 = cos(qJ(3));
t116 = -t128 * g(1) + t130 * g(2);
t118 = t130 * g(1) + t128 * g(2);
t136 = mrSges(4,1) * t116 - mrSges(4,2) * t118 + Ifges(4,3) * qJDD(3);
t103 = m(2) * t117 + qJDD(1) * mrSges(2,1) - t133 * mrSges(2,2) + t137;
t104 = m(2) * t119 - qJDD(1) * mrSges(2,2) + t141 * t133 + t138;
t135 = -t131 * t103 - t129 * t104;
t134 = mrSges(2,1) * t117 + mrSges(3,1) * t113 - mrSges(2,2) * t119 - mrSges(3,3) * t110 + qJ(2) * t137 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);
t132 = qJD(3) ^ 2;
t109 = m(4) * t118 - t132 * mrSges(4,1) - qJDD(3) * mrSges(4,2);
t108 = m(4) * t116 + qJDD(3) * mrSges(4,1) - t132 * mrSges(4,2);
t107 = -mrSges(4,2) * g(3) - mrSges(4,3) * t116 + Ifges(4,5) * qJDD(3) - t132 * Ifges(4,6);
t106 = mrSges(4,1) * g(3) + mrSges(4,3) * t118 + t132 * Ifges(4,5) + Ifges(4,6) * qJDD(3);
t102 = -mrSges(3,2) * t110 - mrSges(2,3) * t117 - t140 * t133 + t139 * qJDD(1) + (-mrSges(2,2) + mrSges(3,1)) * g(3);
t101 = -mrSges(3,2) * t113 + mrSges(2,3) * t119 + t139 * t133 + t140 * qJDD(1) + (m(3) * qJ(2) - t141) * g(3);
t1 = [-m(1) * g(1) + t129 * t103 - t131 * t104 + t128 * t108 - t130 * t109; -m(1) * g(2) - t130 * t108 - t128 * t109 + t135; (-m(1) - m(4) + t142) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t129 * t101 - t131 * t102 + t128 * t106 - t130 * t107; -mrSges(1,3) * g(1) - t131 * t101 - t129 * t102 - t130 * t106 - t128 * t107 + (-pkin(1) * t142 + mrSges(1,1)) * g(3); pkin(1) * t135 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t134 + t136; t134; -t133 * mrSges(3,3) + t138; t136; 0;];
tauJB = t1;
