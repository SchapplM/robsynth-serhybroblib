% Calculate vector of inverse dynamics joint torques for
% fourbar1OL
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:10
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1OL_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:10:25
% EndTime: 2020-04-24 20:10:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (42->34), mult. (74->53), div. (0->0), fcn. (17->6), ass. (0->10)
t11 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1);
t10 = qJD(1) ^ 2;
t8 = cos(qJ(1));
t7 = cos(qJ(2));
t6 = sin(qJ(1));
t5 = sin(qJ(2));
t3 = qJDD(1) + qJDD(2) / 0.2e1;
t2 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t1 = -mrSges(3,1) * g(1) - mrSges(3,2) * g(2);
t4 = [(-g(2) * mrSges(2,1) + mrSges(2,2) * g(1) + t1 * t5 + t2 * t7) * t8 + (g(1) * mrSges(2,1) + mrSges(2,2) * g(2) + t1 * t7 - t2 * t5) * t6 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1) + Ifges(3,3) * qJDD(2) + (0.2e1 * (-t3 * mrSges(3,1) + mrSges(3,2) * t11) * t7 + 0.2e1 * (mrSges(3,1) * t11 + mrSges(3,2) * t3) * t5 + (g(1) * t6 - g(2) * t8 + pkin(2) * qJDD(1)) * m(3)) * pkin(2); (t2 * t8 + t1 * t6 - pkin(2) * (mrSges(3,1) * qJDD(1) + mrSges(3,2) * t10)) * t7 + Ifges(3,3) * (qJDD(1) + qJDD(2)) + (t1 * t8 - t2 * t6 - pkin(2) * (mrSges(3,1) * t10 - mrSges(3,2) * qJDD(1))) * t5; (-g(2) * mrSges(4,1) + g(1) * mrSges(4,2)) * cos(qJ(3)) + (g(1) * mrSges(4,1) + mrSges(4,2) * g(2)) * sin(qJ(3)) + Ifges(4,3) * qJDD(3); 0;];
tau = t4;
