% Calculate vector of inverse dynamics joint torques for with Newton-Euler
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = fourbar2OL_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_invdynJ_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar2OL_invdynJ_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar2OL_invdynJ_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:25
% EndTime: 2020-04-24 20:32:25
% DurationCPUTime: 0.05s
% Computational Cost: add. (42->29), mult. (69->42), div. (0->0), fcn. (36->6), ass. (0->16)
t126 = m(3) * pkin(2);
t117 = sin(qJ(2));
t125 = t117 * mrSges(3,2);
t120 = cos(qJ(2));
t124 = t120 * mrSges(3,1);
t123 = t117 * mrSges(3,1) + t120 * mrSges(3,2);
t122 = -t124 + t125;
t121 = cos(qJ(1));
t119 = cos(qJ(3));
t118 = sin(qJ(1));
t116 = sin(qJ(3));
t112 = -mrSges(2,2) + t123;
t111 = -qJD(1) ^ 2 * pkin(2) - t121 * g(1) - t118 * g(2);
t110 = t118 * g(1) - t121 * g(2) + qJDD(1) * pkin(2);
t109 = mrSges(2,1) + t122 + t126;
t1 = [-(-t118 * t109 + t112 * t121) * g(1) - (t109 * t121 + t118 * t112) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1) + Ifges(3,3) * qJDD(2) + (t122 * qJDD(2) + (-0.2e1 * t124 + 0.2e1 * t125 + t126) * qJDD(1) + (0.2e1 * qJD(1) + qJD(2)) * t123 * qJD(2)) * pkin(2); Ifges(3,3) * (qJDD(1) + qJDD(2)) + mrSges(3,1) * (-t120 * t110 + t117 * t111) - mrSges(3,2) * (-t117 * t110 - t120 * t111); -(-t116 * mrSges(4,1) - mrSges(4,2) * t119) * g(1) - (mrSges(4,1) * t119 - t116 * mrSges(4,2)) * g(2) + Ifges(4,3) * qJDD(3); 0;];
tauJ = t1;
