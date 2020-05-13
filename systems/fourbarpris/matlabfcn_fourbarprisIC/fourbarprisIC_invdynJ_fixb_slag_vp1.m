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
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbarprisIC_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisIC_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:31
% EndTime: 2020-05-07 09:59:31
% DurationCPUTime: 0.37s
% Computational Cost: add. (146->52), mult. (344->74), div. (4->4), fcn. (174->4), ass. (0->24)
t22 = sin(qJ(1));
t38 = t22 * rSges(3,3);
t26 = qJD(1) ^ 2;
t24 = cos(qJ(1));
t36 = t22 * rSges(3,1);
t28 = t24 * rSges(3,3) + t36;
t37 = qJD(1) * t28;
t34 = rSges(3,3) + qJ(2);
t32 = qJD(2) * t22;
t31 = qJ(2) * qJD(1);
t16 = t24 * rSges(3,1) - t38;
t30 = -qJD(2) * t24 + t22 * t31;
t29 = t24 * rSges(2,1) - t22 * rSges(2,2);
t14 = t22 * rSges(2,1) + t24 * rSges(2,2);
t21 = sin(qJ(3));
t23 = cos(qJ(3));
t27 = t23 * rSges(4,1) - t21 * rSges(4,2);
t13 = t21 * rSges(4,1) + t23 * rSges(4,2);
t4 = -qJD(1) * t16 + t30;
t19 = qJD(1) * t38;
t5 = t24 * t31 + t32 + t37;
t2 = -qJDD(2) * t22 - qJDD(1) * t28 + (t19 + (-rSges(3,1) * qJD(1) - 0.2e1 * qJD(2)) * t24) * qJD(1) + (-qJDD(1) * t24 + t26 * t22) * qJ(2);
t1 = -qJDD(1) * t16 - qJDD(2) * t24 + (qJDD(1) * t22 + t26 * t24) * qJ(2) + (0.2e1 * t32 + t37) * qJD(1);
t3 = [(-t23 * t21 - t24 * t22) / (pkin(3) + qJ(2)) / (-t23 ^ 2 + t24 ^ 2) * ((Icges(2,3) + Icges(3,2)) * qJDD(1) + (t4 * t32 - t5 * (t19 + t30) + (t4 * t36 + (t5 * rSges(3,1) + t4 * t34) * t24) * qJD(1) + (t2 - g(2)) * (-t34 * t24 - t36) + (t1 - g(1)) * (t22 * qJ(2) - t16)) * m(3) + ((qJDD(1) * t29 + g(2)) * t29 + (qJDD(1) * t14 - g(1)) * t14) * m(2)) + m(3) * (-t1 * t24 - t2 * t22) - m(3) * (-g(1) * t24 - g(2) * t22) - 0.1e1 / pkin(2) / (t21 * t24 - t22 * t23) * (Icges(4,3) * qJDD(3) + (-(-qJDD(3) * t27 - g(2)) * t27 + (qJDD(3) * t13 - g(1)) * t13) * m(4));];
tau = t3(:);
