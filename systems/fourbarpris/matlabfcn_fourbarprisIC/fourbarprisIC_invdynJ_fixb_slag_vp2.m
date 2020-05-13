% Calculate vector of inverse dynamics joint torques with ic for
% fourbarprisIC
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
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbarprisIC_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:32
% EndTime: 2020-05-07 09:59:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (34->27), mult. (58->36), div. (4->4), fcn. (26->4), ass. (0->9)
t3 = sin(qJ(1));
t5 = cos(qJ(1));
t13 = g(1) * t5 + g(2) * t3 + qJDD(2);
t10 = m(3) * qJ(2);
t9 = qJ(2) * qJDD(1);
t7 = -g(1) * t3 + g(2) * t5 + (2 * qJD(1) * qJD(2)) + t9;
t4 = cos(qJ(3));
t2 = sin(qJ(3));
t1 = [(-t4 * t2 - t5 * t3) / (pkin(3) + qJ(2)) / (-t4 ^ 2 + t5 ^ 2) * (-g(1) * (t3 * mrSges(2,1) + t5 * mrSges(2,2)) - g(2) * (-t5 * mrSges(2,1) + t3 * mrSges(2,2)) + (Ifges(2,3) + Ifges(3,2)) * qJDD(1) + t13 * mrSges(3,1) + t7 * t10 + (t7 + t9) * mrSges(3,3)) + qJDD(1) * mrSges(3,1) - 0.1e1 / pkin(2) / (t2 * t5 - t3 * t4) * (Ifges(4,3) * qJDD(3) - g(1) * (t2 * mrSges(4,1) + t4 * mrSges(4,2)) - g(2) * (-t4 * mrSges(4,1) + t2 * mrSges(4,2))) + (-mrSges(3,3) - t10) * (qJD(1) ^ 2) + t13 * m(3);];
tau = t1(:);
