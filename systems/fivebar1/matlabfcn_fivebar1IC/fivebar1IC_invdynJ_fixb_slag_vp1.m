% Calculate vector of inverse dynamics joint torques with ic for
% fivebar1IC
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
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fivebar1IC_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:09
% EndTime: 2020-04-27 06:19:10
% DurationCPUTime: 0.27s
% Computational Cost: add. (140->79), mult. (296->118), div. (8->4), fcn. (44->19), ass. (0->30)
t44 = (rSges(5,1) ^ 2 + rSges(5,2) ^ 2) * m(5);
t43 = (rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * m(3);
t42 = pkin(2) * m(3);
t41 = pkin(3) * m(5);
t19 = qJ(3) + qJ(4);
t40 = ((rSges(5,1) * g(1) + rSges(5,2) * g(2)) * sin(t19) + (-rSges(5,1) * g(2) + rSges(5,2) * g(1)) * cos(t19)) * m(5);
t20 = qJ(1) + qJ(2);
t39 = ((-rSges(3,1) * g(1) - rSges(3,2) * g(2)) * sin(t20) + (rSges(3,1) * g(2) - rSges(3,2) * g(1)) * cos(t20)) * m(3);
t17 = qJDD(3) + qJDD(4);
t21 = sin(qJ(4));
t23 = cos(qJ(4));
t25 = qJD(3) ^ 2;
t38 = ((Icges(5,3) + t44) * t17 + ((rSges(5,1) * qJDD(3) + rSges(5,2) * t25) * t23 + (rSges(5,1) * t25 - rSges(5,2) * qJDD(3)) * t21) * t41 + t40) / pkin(4);
t18 = qJDD(1) + qJDD(2);
t22 = sin(qJ(2));
t24 = cos(qJ(2));
t26 = qJD(1) ^ 2;
t37 = (t18 * (Icges(3,3) + t43) + (-(rSges(3,1) * qJDD(1) + rSges(3,2) * t26) * t24 - (rSges(3,1) * t26 - rSges(3,2) * qJDD(1)) * t22) * t42 + t39) / pkin(5);
t36 = -cos(0.2e1 * t19) + cos(0.2e1 * t20);
t35 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1);
t34 = qJD(4) * (qJD(3) + qJD(4) / 0.2e1);
t33 = 0.2e1 * t42;
t32 = -0.2e1 * t41;
t28 = 2 * qJ(1);
t27 = 2 * qJ(3);
t12 = qJDD(1) + qJDD(2) / 0.2e1;
t11 = qJDD(3) + qJDD(4) / 0.2e1;
t8 = -0.1e1 / sin(t19 - t20);
t3 = 0.1e1 / t36;
t1 = [(-t12 * rSges(3,1) + rSges(3,2) * t35) * t24 * t33 + t18 * t43 + Icges(3,3) * qJDD(2) + (-t36 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t27) + cos(qJ(2) + t28)) * pkin(2)) * t3 * t37 + t39 + (-g(2) * t42 + (-rSges(2,1) * g(2) + rSges(2,2) * g(1)) * m(2)) * cos(qJ(1)) + (g(1) * t42 + (rSges(2,1) * g(1) + rSges(2,2) * g(2)) * m(2)) * sin(qJ(1)) + ((rSges(3,1) * t35 + rSges(3,2) * t12) * t33 + pkin(2) * t8 * t38) * t22 + (pkin(2) ^ 2 * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3)) * qJDD(1); (-t11 * rSges(5,1) + rSges(5,2) * t34) * t23 * t32 + t17 * t44 + Icges(5,3) * qJDD(4) + (-t36 * pkin(4) + (cos(qJ(4) + t27) - cos(qJ(4) + t28 + 0.2e1 * qJ(2))) * pkin(3)) * t3 * t38 + t40 + (-g(2) * t41 + (-rSges(4,1) * g(2) + rSges(4,2) * g(1)) * m(4)) * cos(qJ(3)) + (g(1) * t41 + (rSges(4,1) * g(1) + rSges(4,2) * g(2)) * m(4)) * sin(qJ(3)) + (pkin(3) * t8 * t37 + (rSges(5,1) * t34 + rSges(5,2) * t11) * t32) * t21 + (pkin(3) ^ 2 * m(5) + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3) + Icges(5,3)) * qJDD(3);];
tau = t1(:);
