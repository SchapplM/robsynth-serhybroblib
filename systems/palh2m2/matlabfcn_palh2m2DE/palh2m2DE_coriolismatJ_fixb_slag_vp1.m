% Calculate matrix of centrifugal and coriolis load on the joints for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh2m2DE_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:28
% EndTime: 2020-05-03 01:06:29
% DurationCPUTime: 0.32s
% Computational Cost: add. (283->76), mult. (773->113), div. (0->0), fcn. (209->6), ass. (0->55)
t30 = sin(qJ(2));
t74 = pkin(4) * t30;
t29 = sin(qJ(3));
t73 = pkin(5) * t29;
t28 = sin(qJ(4));
t59 = rSges(7,2) * t28;
t31 = cos(qJ(4));
t60 = rSges(7,1) * t31;
t25 = -t59 + t60;
t39 = t25 * m(7);
t72 = t39 / 0.2e1 + m(7) * (-t60 / 0.2e1 + t59 / 0.2e1);
t48 = pkin(1) + pkin(2);
t24 = rSges(7,1) * t28 + rSges(7,2) * t31;
t27 = pkin(3) + t48;
t33 = cos(qJ(2));
t61 = t33 * pkin(4);
t32 = cos(qJ(3));
t62 = t32 * pkin(5);
t71 = (t27 + t61 + t62) * t24;
t70 = m(5) + m(6);
t69 = m(6) + m(7);
t68 = m(5) * rSges(5,1);
t67 = m(5) * rSges(5,2);
t66 = m(6) * rSges(6,1);
t65 = m(7) * t24;
t64 = rSges(3,3) * m(3);
t63 = rSges(5,3) * m(5);
t26 = pkin(5) * t69 + t68;
t46 = t32 * t67;
t58 = (t26 * t29 + t46) * t33;
t57 = pkin(4) * qJD(2);
t56 = m(3) * pkin(1);
t55 = m(7) * qJD(4);
t40 = t25 + t27;
t42 = pkin(5) * m(6) + t68;
t43 = -rSges(5,1) * t67 + Icges(5,4);
t1 = t48 * t46 + pkin(4) * t58 + (t43 * t29 + (m(7) * t40 + t66) * pkin(5) + t48 * t42) * t29 + (((rSges(5,1) ^ 2 - rSges(5,2) ^ 2) * m(5) - Icges(5,1) + Icges(5,2) + t69 * pkin(5) ^ 2) * t29 - t43 * t32) * t32;
t54 = t1 * qJD(1);
t44 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t47 = t29 * t67;
t49 = m(4) + t70;
t2 = (rSges(3,2) * t56 - t44 * t33) * t33 + (((rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) + Icges(3,2) + (m(7) + t49) * pkin(4) ^ 2) * t33 + rSges(3,1) * t56 + t44 * t30 + (t42 * t32 - t47 + t49 * pkin(1) + t70 * pkin(2) + m(4) * rSges(4,1) + t66 + (t40 + t62) * m(7)) * pkin(4)) * t30;
t53 = t2 * qJD(1);
t5 = m(7) * t71;
t52 = t5 * qJD(1);
t6 = t39 * t74;
t51 = t6 * qJD(1);
t7 = t39 * t73;
t50 = t7 * qJD(1);
t45 = t25 * t55;
t41 = -(t26 * t32 - t47) * t30 + t58;
t36 = -rSges(6,3) * m(6) + t65;
t9 = t72 * t74;
t8 = t72 * t73;
t3 = [-qJD(2) * t2 - qJD(3) * t1 - qJD(4) * t5, -t53 + ((-rSges(3,1) * t64 + Icges(3,5)) * t33 + (rSges(3,2) * t64 - Icges(3,6)) * t30 + (-rSges(4,3) * m(4) + t36 - t63) * t61) * qJD(2) + t9 * qJD(4), -t54 + ((-rSges(5,1) * t63 + Icges(5,5)) * t32 - (-rSges(5,2) * t63 + Icges(5,6)) * t29 + t36 * t62) * qJD(3) + t8 * qJD(4), t9 * qJD(2) + t8 * qJD(3) - t55 * t71 - t52; qJD(4) * t6 + t53, 0, -t41 * pkin(4) * qJD(3), t45 * t74 + t51; qJD(4) * t7 + t54, t41 * t57, 0, t45 * t73 + t50; -qJD(2) * t6 - qJD(3) * t7 + t52, t33 * t57 * t65 - t51, qJD(3) * t62 * t65 - t50, 0;];
Cq = t3;
