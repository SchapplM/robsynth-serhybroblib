% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
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
% Cq [1x1]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1TE_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1TE_coriolismatJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:49:17
% EndTime: 2020-04-24 19:49:19
% DurationCPUTime: 0.70s
% Computational Cost: add. (1889->73), mult. (2622->148), div. (84->11), fcn. (691->4), ass. (0->69)
t92 = pkin(4) ^ 2;
t54 = pkin(2) ^ 2;
t55 = pkin(1) ^ 2;
t48 = cos(qJ(1));
t81 = pkin(2) * t48;
t73 = -0.2e1 * pkin(1) * t81 + t55;
t38 = t54 + t73;
t72 = pkin(3) ^ 2 - t92;
t27 = t38 - t72;
t39 = pkin(1) - t81;
t86 = -pkin(3) - pkin(4);
t20 = (pkin(2) - t86) * (pkin(2) + t86) + t73;
t85 = pkin(4) - pkin(3);
t21 = (pkin(2) - t85) * (pkin(2) + t85) + t73;
t75 = t20 * t21;
t56 = sqrt(-t75);
t47 = sin(qJ(1));
t80 = t47 * pkin(2);
t11 = t27 * t80 + t39 * t56;
t15 = 0.1e1 / t56;
t32 = 0.1e1 / t38;
t83 = pkin(1) * t32;
t70 = t15 * t83;
t3 = -t11 * t70 + 0.1e1;
t52 = 0.1e1 / pkin(3);
t91 = t3 * t52;
t18 = 0.1e1 / t21;
t90 = t18 / t20 ^ 2;
t58 = t38 ^ 2;
t33 = 0.1e1 / t58;
t34 = 0.1e1 / t58 ^ 2;
t88 = -t32 / 0.2e1;
t87 = t32 / 0.2e1;
t78 = rSges(4,2) * t47;
t22 = t39 * rSges(4,1) + pkin(2) * t78;
t31 = rSges(4,1) * t80 + rSges(4,2) * t81;
t23 = -rSges(4,2) * pkin(1) + t31;
t6 = t22 * t56 + t27 * t23;
t7 = -t27 * t22 + t23 * t56;
t84 = t6 ^ 2 + t7 ^ 2;
t82 = pkin(1) * t47;
t79 = rSges(3,2) * t47;
t71 = pkin(1) * t80;
t13 = (-t20 - t21) * t71;
t77 = t13 * t15;
t16 = 0.1e1 / t20;
t76 = t16 * t18;
t74 = t47 * t56;
t30 = rSges(3,1) * t80 + rSges(3,2) * t81;
t69 = pkin(2) * t54 * t82;
t24 = -rSges(3,2) * pkin(1) + t30;
t25 = t39 * rSges(3,1) + pkin(2) * t79;
t26 = t38 + t72;
t8 = t24 * t56 + t26 * t25;
t68 = t8 * t87;
t9 = -t26 * t24 + t25 * t56;
t67 = t9 * t88;
t40 = pkin(1) * t48 - pkin(2);
t12 = -t26 * t82 + t40 * t56;
t46 = t47 ^ 2;
t66 = t12 * (t40 * t77 - 0.2e1 * t55 * t46 * pkin(2) + (-t48 * t26 - t74) * pkin(1)) * t54 * t76;
t65 = t84 * t69;
t63 = t34 * t65;
t29 = (rSges(3,1) * t48 - t79) * pkin(2);
t28 = (rSges(4,1) * t48 - t78) * pkin(2);
t19 = 0.1e1 / t21 ^ 2;
t10 = t12 ^ 2;
t1 = -(t39 * t77 + 0.2e1 * t54 * t46 * pkin(1) + (t48 * t27 + t74) * pkin(2)) * t70 + (-t13 / t75 * t83 + 0.2e1 * t55 * t33 * t80) * t11 * t15;
t2 = [(t3 * Icges(3,3) * t1 + m(3) * ((t67 * t91 - t80) * (-t81 + (t1 * t67 + ((t25 * t77 - t26 * t29 + t30 * t56) * t88 + (t24 * t32 + t33 * t9) * t71) * t3) * t52) + (t68 * t91 + t81) * (-t80 + (t1 * t68 + ((t24 * t77 + t26 * t30 + t29 * t56) * t87 + (t25 * t32 - t33 * t8) * t71) * t3) * t52)) + m(4) * (-t84 * t34 * t66 + (t63 * t90 + (t19 * t63 + (0.4e1 * t32 * t65 + (-t6 * (t22 * t77 + 0.2e1 * t23 * t71 + t27 * t28 + t31 * t56) - t7 * (-0.2e1 * t22 * t71 + t23 * t77 - t27 * t31 + t28 * t56)) * t54) * t18 * t34) * t16) * t10) / t92 / 0.4e1 + (-t66 + (t16 * t19 + 0.2e1 * t32 * t76 + t90) * t10 * t69) * t33 * Icges(4,3)) * qJD(1);];
Cq = t2;
