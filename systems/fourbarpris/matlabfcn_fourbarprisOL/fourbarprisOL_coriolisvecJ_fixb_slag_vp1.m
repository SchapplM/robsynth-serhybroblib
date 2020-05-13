% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbarprisOL_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_coriolisvecJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisOL_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisOL_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:05
% EndTime: 2020-05-07 09:52:06
% DurationCPUTime: 0.13s
% Computational Cost: add. (91->24), mult. (254->38), div. (0->0), fcn. (116->2), ass. (0->16)
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t33 = (t18 * rSges(3,1) + t20 * rSges(3,3)) * qJD(1);
t16 = t18 * rSges(3,3);
t32 = qJD(1) ^ 2 * qJ(2);
t31 = rSges(3,3) + qJ(2);
t30 = rSges(3,1) * qJD(1);
t29 = qJD(2) * t18;
t28 = qJ(2) * qJD(1);
t14 = t18 * t28;
t4 = -qJD(1) * (t20 * rSges(3,1) - t16) - qJD(2) * t20 + t14;
t15 = qJD(1) * t16;
t5 = t20 * t28 + t29 + t33;
t3 = t20 * t32 + (0.2e1 * t29 + t33) * qJD(1);
t2 = t18 * t32 + (t15 + (-0.2e1 * qJD(2) - t30) * t20) * qJD(1);
t1 = [m(3) * (t3 * t16 - t5 * (t14 + t15) + (t3 * qJ(2) + t4 * (qJD(2) + t30) - t2 * rSges(3,1)) * t18 + (-t3 * rSges(3,1) - t2 * t31 + t5 * qJD(2) + (t5 * rSges(3,1) + t4 * t31) * qJD(1)) * t20); m(3) * (-t2 * t18 - t3 * t20); 0; 0;];
tauc = t1(:);
