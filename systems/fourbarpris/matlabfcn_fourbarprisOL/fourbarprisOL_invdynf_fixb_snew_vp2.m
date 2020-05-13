% Calculate vector of cutting forces with Newton-Euler
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
% f_new [3x4]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = fourbarprisOL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisOL_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_invdynf_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_invdynf_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisOL_invdynf_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisOL_invdynf_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:05
% EndTime: 2020-05-07 09:52:06
% DurationCPUTime: 0.07s
% Computational Cost: add. (95->36), mult. (158->46), div. (0->0), fcn. (48->4), ass. (0->17)
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t16 = -t11 * g(1) + t13 * g(2);
t21 = m(3) * (-qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t16);
t20 = -m(2) - m(3);
t15 = qJD(1) ^ 2;
t18 = t13 * g(1) + t11 * g(2);
t19 = m(3) * (-t15 * qJ(2) + qJDD(2) + t18) + qJDD(1) * mrSges(3,1);
t17 = (mrSges(2,1) + mrSges(3,3));
t14 = qJD(3) ^ 2;
t12 = cos(qJ(3));
t10 = sin(qJ(3));
t4 = m(4) * (t12 * g(1) + t10 * g(2)) - qJDD(3) * mrSges(4,2) - t14 * mrSges(4,1);
t3 = m(4) * (-t10 * g(1) + t12 * g(2)) + qJDD(3) * mrSges(4,1) - t14 * mrSges(4,2);
t2 = m(2) * t18 - qJDD(1) * mrSges(2,2) - (t17 * t15) + t19;
t1 = m(2) * t16 - t21 + (-mrSges(2,2) + mrSges(3,1)) * t15 + t17 * qJDD(1);
t5 = [-m(1) * g(1) + t11 * t1 + t10 * t3 - t12 * t4 - t13 * t2, t2, -t15 * mrSges(3,1) - qJDD(1) * mrSges(3,3) + t21, t4, 0, 0; -m(1) * g(2) - t13 * t1 - t10 * t4 - t11 * t2 - t12 * t3, t1, m(3) * g(3), t3, 0, 0; (-m(1) - m(4) + t20) * g(3), t20 * g(3), -t15 * mrSges(3,3) + t19, -m(4) * g(3), 0, 0;];
f_new = t5;
