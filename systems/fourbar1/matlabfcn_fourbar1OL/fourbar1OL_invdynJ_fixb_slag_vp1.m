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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:10
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1OL_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1OL_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:10:25
% EndTime: 2020-04-24 20:10:25
% DurationCPUTime: 0.11s
% Computational Cost: add. (47->37), mult. (114->61), div. (0->0), fcn. (12->8), ass. (0->11)
t16 = pkin(2) * m(3);
t6 = qJDD(1) + qJDD(2);
t15 = (rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * t6;
t7 = qJ(1) + qJ(2);
t14 = ((-g(1) * rSges(3,1) - rSges(3,2) * g(2)) * sin(t7) + (g(2) * rSges(3,1) - rSges(3,2) * g(1)) * cos(t7)) * m(3);
t13 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1);
t10 = qJD(1) ^ 2;
t9 = cos(qJ(2));
t8 = sin(qJ(2));
t3 = qJDD(1) + qJDD(2) / 0.2e1;
t1 = [t15 * m(3) + Icges(3,3) * qJDD(2) + t14 + 0.2e1 * ((-t3 * rSges(3,1) + rSges(3,2) * t13) * t9 + (rSges(3,1) * t13 + rSges(3,2) * t3) * t8) * t16 + (-g(2) * t16 + m(2) * (-g(2) * rSges(2,1) + rSges(2,2) * g(1))) * cos(qJ(1)) + (g(1) * t16 + m(2) * (g(1) * rSges(2,1) + rSges(2,2) * g(2))) * sin(qJ(1)) + (pkin(2) ^ 2 * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3)) * qJDD(1); Icges(3,3) * t6 + (t15 + (-(rSges(3,1) * qJDD(1) + rSges(3,2) * t10) * t9 - (rSges(3,1) * t10 - rSges(3,2) * qJDD(1)) * t8) * pkin(2)) * m(3) + t14; Icges(4,3) * qJDD(3) + ((-g(2) * rSges(4,1) + rSges(4,2) * g(1)) * cos(qJ(3)) + (rSges(4,1) * g(1) + rSges(4,2) * g(2)) * sin(qJ(3)) + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * qJDD(3)) * m(4); 0;];
tau = t1;
