% Calculate vector of cutting forces with Newton-Euler
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
% f_new [3x4]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:10
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = fourbar1OL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1OL_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar1OL_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1OL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_invdynf_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1OL_invdynf_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1OL_invdynf_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1OL_invdynf_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:10:27
% EndTime: 2020-04-24 20:10:27
% DurationCPUTime: 0.22s
% Computational Cost: add. (130->56), mult. (196->76), div. (0->0), fcn. (108->6), ass. (0->29)
t34 = (2 * qJD(1) + qJD(2)) * qJD(2);
t26 = m(2) + m(3);
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t9 = t21 * mrSges(3,1) + t24 * mrSges(3,2);
t10 = t24 * mrSges(3,1) - t21 * mrSges(3,2);
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t32 = -t25 * g(1) - t22 * g(2);
t31 = t22 * g(1) - t25 * g(2);
t30 = t10 * t25 - t22 * t9;
t4 = t22 * t10 + t9 * t25;
t29 = qJD(1) ^ 2;
t27 = qJD(3) ^ 2;
t23 = cos(qJ(3));
t20 = sin(qJ(3));
t19 = qJD(1) + qJD(2);
t18 = qJDD(1) + qJDD(2);
t17 = t19 ^ 2;
t16 = m(4) + m(1) + t26;
t12 = mrSges(4,1) * t23 - t20 * mrSges(4,2);
t11 = -t20 * mrSges(4,1) - mrSges(4,2) * t23;
t8 = -mrSges(2,2) + t9;
t7 = -(t29 * pkin(2)) + t32;
t6 = (qJDD(1) * pkin(2)) + t31;
t5 = -(m(3) * pkin(2)) - mrSges(2,1) + t10;
t2 = t22 * t5 + t8 * t25;
t1 = t22 * t8 - t5 * t25;
t3 = [-t16 * g(1) + t2 * qJDD(1) + t4 * qJDD(2) + t11 * qJDD(3) - t1 * t29 - t12 * t27 + t34 * t30, t8 * qJDD(1) + t9 * qJDD(2) + t34 * t10 + t32 * t26 + t5 * t29, m(3) * (-t21 * t6 - t24 * t7) - t18 * mrSges(3,2) - t17 * mrSges(3,1), -mrSges(4,1) * t27 - mrSges(4,2) * qJDD(3) + (-g(1) * t23 - g(2) * t20) * m(4), 0, 0; -t16 * g(2) + t1 * qJDD(1) - t30 * qJDD(2) + t12 * qJDD(3) + t11 * t27 + t2 * t29 + t34 * t4, -t5 * qJDD(1) - t10 * qJDD(2) + t31 * t26 + t8 * t29 + t34 * t9, m(3) * (t21 * t7 - t24 * t6) + t18 * mrSges(3,1) - t17 * mrSges(3,2), mrSges(4,1) * qJDD(3) - mrSges(4,2) * t27 + (g(1) * t20 - g(2) * t23) * m(4), 0, 0; -t16 * g(3), -t26 * g(3), -m(3) * g(3), -m(4) * g(3), 0, 0;];
f_new = t3;
