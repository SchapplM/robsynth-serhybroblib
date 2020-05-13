% Calculate vector of inverse dynamics joint torques for
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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fivebar1OL_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:09
% EndTime: 2020-04-27 06:13:10
% DurationCPUTime: 0.18s
% Computational Cost: add. (82->62), mult. (200->98), div. (0->0), fcn. (20->12), ass. (0->21)
t31 = (rSges(5,1) ^ 2 + rSges(5,2) ^ 2) * m(5);
t30 = (rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * m(3);
t29 = pkin(2) * m(3);
t28 = pkin(3) * m(5);
t13 = qJ(3) + qJ(4);
t27 = ((g(1) * rSges(5,1) + rSges(5,2) * g(2)) * sin(t13) + (-g(2) * rSges(5,1) + rSges(5,2) * g(1)) * cos(t13)) * m(5);
t14 = qJ(1) + qJ(2);
t26 = ((-g(1) * rSges(3,1) - rSges(3,2) * g(2)) * sin(t14) + (g(2) * rSges(3,1) - rSges(3,2) * g(1)) * cos(t14)) * m(3);
t25 = qJD(4) * (qJD(3) + qJD(4) / 0.2e1);
t24 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1);
t20 = qJD(1) ^ 2;
t19 = qJD(3) ^ 2;
t18 = cos(qJ(2));
t17 = cos(qJ(4));
t16 = sin(qJ(2));
t15 = sin(qJ(4));
t12 = qJDD(1) + qJDD(2);
t11 = qJDD(3) + qJDD(4);
t6 = qJDD(1) + qJDD(2) / 0.2e1;
t5 = qJDD(3) + qJDD(4) / 0.2e1;
t1 = [t12 * t30 + Icges(3,3) * qJDD(2) + t26 + 0.2e1 * ((-t6 * rSges(3,1) + rSges(3,2) * t24) * t18 + (rSges(3,1) * t24 + rSges(3,2) * t6) * t16) * t29 + (-g(2) * t29 + m(2) * (-g(2) * rSges(2,1) + rSges(2,2) * g(1))) * cos(qJ(1)) + (g(1) * t29 + m(2) * (g(1) * rSges(2,1) + rSges(2,2) * g(2))) * sin(qJ(1)) + (m(3) * pkin(2) ^ 2 + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3)) * qJDD(1); t12 * (Icges(3,3) + t30) + (-(rSges(3,1) * qJDD(1) + rSges(3,2) * t20) * t18 - (rSges(3,1) * t20 - rSges(3,2) * qJDD(1)) * t16) * t29 + t26; t11 * t31 + Icges(5,3) * qJDD(4) + t27 - 0.2e1 * ((-t5 * rSges(5,1) + rSges(5,2) * t25) * t17 + (rSges(5,1) * t25 + rSges(5,2) * t5) * t15) * t28 + (-g(2) * t28 + m(4) * (-g(2) * rSges(4,1) + rSges(4,2) * g(1))) * cos(qJ(3)) + (g(1) * t28 + m(4) * (g(1) * rSges(4,1) + rSges(4,2) * g(2))) * sin(qJ(3)) + (pkin(3) ^ 2 * m(5) + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3) + Icges(5,3)) * qJDD(3); (Icges(5,3) + t31) * t11 + ((rSges(5,1) * qJDD(3) + rSges(5,2) * t19) * t17 + (rSges(5,1) * t19 - rSges(5,2) * qJDD(3)) * t15) * t28 + t27; 0;];
tau = t1;
