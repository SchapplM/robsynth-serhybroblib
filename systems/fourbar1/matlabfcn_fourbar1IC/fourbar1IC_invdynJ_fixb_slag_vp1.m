% Calculate vector of inverse dynamics joint torques with ic for
% fourbar1IC
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:15
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1IC_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1IC_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1IC_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1IC_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1IC_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:15:49
% EndTime: 2020-04-24 20:15:49
% DurationCPUTime: 0.20s
% Computational Cost: add. (57->42), mult. (119->66), div. (4->3), fcn. (17->10), ass. (0->13)
t19 = pkin(2) * m(3);
t8 = qJDD(1) + qJDD(2);
t18 = (rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * t8;
t9 = qJ(1) + qJ(2);
t17 = ((-rSges(3,1) * g(1) - rSges(3,2) * g(2)) * sin(t9) + (rSges(3,1) * g(2) - rSges(3,2) * g(1)) * cos(t9)) * m(3);
t16 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1);
t14 = -qJ(3) + qJ(1);
t12 = qJD(1) ^ 2;
t11 = cos(qJ(2));
t10 = sin(qJ(2));
t5 = qJDD(1) + qJDD(2) / 0.2e1;
t4 = sin(qJ(2) + t14);
t1 = [t18 * m(3) + Icges(3,3) * qJDD(2) + t17 + (pkin(2) ^ 2 * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3)) * qJDD(1) + ((-pkin(3) * t4 + pkin(2) * sin(t14)) / pkin(3) * (Icges(3,3) * t8 + (t18 + (-(rSges(3,1) * qJDD(1) + rSges(3,2) * t12) * t11 - (rSges(3,1) * t12 - rSges(3,2) * qJDD(1)) * t10) * pkin(2)) * m(3) + t17) + t10 * pkin(2) / pkin(4) * (Icges(4,3) * qJDD(3) + ((-rSges(4,1) * g(2) + rSges(4,2) * g(1)) * cos(qJ(3)) + (rSges(4,1) * g(1) + rSges(4,2) * g(2)) * sin(qJ(3)) + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * qJDD(3)) * m(4))) / t4 + 0.2e1 * ((-t5 * rSges(3,1) + rSges(3,2) * t16) * t11 + (rSges(3,1) * t16 + rSges(3,2) * t5) * t10) * t19 + (-g(2) * t19 + (-rSges(2,1) * g(2) + rSges(2,2) * g(1)) * m(2)) * cos(qJ(1)) + (g(1) * t19 + (rSges(2,1) * g(1) + rSges(2,2) * g(2)) * m(2)) * sin(qJ(1));];
tau = t1(:);
